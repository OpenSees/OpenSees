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
                                                                        
// $Revision: 1.33 $
// $Date: 2003-05-15 21:30:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/element/NLBeamColumn2d.cpp,v $
                                                                        
                                                                        
//
// Written by Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu) on 01/99 
// Revised: rms 01/99 (distributed loads)
//          rms 06/99 (mass matrix)
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
//          rms 07/00 (new state determination procedure, no need to store fscommit)
//
// Purpose: This file contains the implementation for the NLBeamColumn2d class.
//          NLBeamColumn2d.C is a materially nonlinear flexibility based frame element.

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <Information.h>
#include <NLBeamColumn2d.h>
#include <MatrixUtil.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <math.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>

#define  NDM   2         // dimension of the problem (2d)
#define  NND   3         // number of nodal dof's
#define  NEGD  6         // number of element global dof's
#define  NEBD  3         // number of element dof's in the basic system

Matrix NLBeamColumn2d::theMatrix(6,6);
Vector NLBeamColumn2d::theVector(6);
double NLBeamColumn2d::workArea[100];
GaussLobattoQuadRule1d01 NLBeamColumn2d::quadRule;
Vector *NLBeamColumn2d::vsSubdivide  = 0;
Matrix *NLBeamColumn2d::fsSubdivide  = 0;
Vector *NLBeamColumn2d::SsrSubdivide = 0;
int NLBeamColumn2d::maxNumSections   = 0;


// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
NLBeamColumn2d::NLBeamColumn2d(): 
Element(0,ELE_TAG_NLBeamColumn2d), connectedExternalNodes(2), 
nSections(0), sections(0), crdTransf(0),
rho(0.0), maxIters(0), tol(0.0),
cosTheta(0.0), sinTheta(0.0), initialFlag(0),
load(NEGD), 
kv(NEBD,NEBD), Se(NEBD),
kvcommit(NEBD,NEBD), Secommit(NEBD),
fs(0), vs(0), Ssr(0), vscommit(0), sp(0), Ki(0), maxSubdivisions(0)
{
  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  
  for (int i=0; i<2; i++)
    theNodes[i] = 0;

}


// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points. 
// allocates the necessary space needed by each object
NLBeamColumn2d::NLBeamColumn2d (int tag, int nodeI, int nodeJ,
                                int numSections, SectionForceDeformation *sectionPtrs[],
                                CrdTransf2d &coordTransf, double massDensPerUnitLength,
				int maxNumIters, double tolerance, int maxSub):
Element(tag,ELE_TAG_NLBeamColumn2d), connectedExternalNodes(2),
nSections(numSections), sections(sectionPtrs), crdTransf(0),
rho(massDensPerUnitLength),maxIters(maxNumIters), tol(tolerance), 
cosTheta(0.0), sinTheta(1.0), initialFlag(0),
load(NEGD), 
kv(NEBD,NEBD), Se(NEBD), 
kvcommit(NEBD,NEBD), Secommit(NEBD),
fs(0), vs(0),Ssr(0), vscommit(0), sp(0), Ki(0), maxSubdivisions(maxSub)
{
   connectedExternalNodes(0) = nodeI;
   connectedExternalNodes(1) = nodeJ;    

   // get copy of the sections
   
   if (!sectionPtrs)
   {
       opserr << "Error: NLBeamColumn2d::NLBeamColumn2d:  invalid section pointer ";
       exit(-1);
   }	  
   
   sections = new SectionForceDeformation *[nSections];
   if (!sections)
   {
       opserr << "Error: NLBeamColumn2d::NLBeamColumn2d: could not alocate section pointer";
       exit(-1);
   }  
   
   for (int i = 0; i < nSections; i++)
   {
      if (!sectionPtrs[i])
      {
	  opserr << "Error: NLBeamColumn2d::NLBeamColumn2d: section pointer " << i << endln;
          exit(-1);
      }  
       
      sections[i] = sectionPtrs[i]->getCopy();
      if (!sections[i])
      {
	  opserr << "Error: NLBeamColumn2d::NLBeamColumn2d: could not create copy of section " << i << endln;
          exit(-1);
      }
   }

   // get copy of the transformation object   
   
   crdTransf = coordTransf.getCopy(); 
   if (!crdTransf)
   {
      opserr << "Error: NLBeamColumn2d::NLBeamColumn2d: could not create copy of coordinate transformation object" << endln;
      exit(-1);
   }

   // alocate section flexibility matrices and section deformation vectors
   fs  = new Matrix [nSections];
   if (!fs)
   {
       opserr << "NLBeamColumn2d::NLBeamColumn2d() -- failed to allocate fs array";
       exit(-1);
   }
   
   vs = new Vector [nSections];
   if (!vs)
   {
       opserr << "NLBeamColumn2d::NLBeamColumn2d() -- failed to allocate vs array";
       exit(-1);
   }

   Ssr  = new Vector [nSections];
   if (!Ssr)
   {
       opserr << "NLBeamColumn2d::NLBeamColumn2d() -- failed to allocate Ssr array";
       exit(-1);
   }
   
   vscommit = new Vector [nSections];
   if (!vscommit)
   {
       opserr << "NLBeamColumn2d::NLBeamColumn2d() -- failed to allocate vscommit array";   
       exit(-1);
   }

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  if (nSections > maxNumSections) {
    maxNumSections = nSections;
    if (vsSubdivide != 0)
      delete [] vsSubdivide;
    vsSubdivide  = new Vector [nSections];
    if (fsSubdivide != 0)
      delete [] fsSubdivide;
    fsSubdivide  = new Matrix [nSections];
    if (SsrSubdivide != 0)
      delete [] SsrSubdivide;
    SsrSubdivide  = new Vector [nSections];
    if (!vsSubdivide || !fsSubdivide || !SsrSubdivide) {
       opserr << "NLBeamColumn2d::NLBeamColumn2d() -- failed to allocate Subdivide arrays";   
       exit(-1);
   }
  }
}



// ~NLBeamColumn2d():
// 	destructor
//      delete must be invoked on any objects created by the object
NLBeamColumn2d::~NLBeamColumn2d()
{
   if (sections) {
      for (int i=0; i < nSections; i++)
         if (sections[i])
            delete sections[i];
      delete [] sections;
   }

   if (fs) 
     delete [] fs;

   if (vs) 
     delete [] vs;

   if (Ssr) 
     delete [] Ssr;

   if (vscommit) 
     delete [] vscommit;

   if (crdTransf)
     delete crdTransf;   

   if (sp != 0)
     delete sp;
   
   if (Ki != 0)
     delete Ki;
}


int
NLBeamColumn2d::getNumExternalNodes(void) const
{
  return 2;
}


const ID &
NLBeamColumn2d::getExternalNodes(void) 
{
   return connectedExternalNodes;
}


Node **
NLBeamColumn2d::getNodePtrs(void) 
{
  return theNodes;
}

int
NLBeamColumn2d::getNumDOF(void) 
{
   return NEGD;
}



void
NLBeamColumn2d::setDomain(Domain *theDomain)
{
   // check Domain is not null - invoked when object removed from a domain
   if (theDomain == 0)
   {
      theNodes[0] = 0;
      theNodes[1] = 0;
        
      opserr << "NLBeamColumn2d::setDomain:  theDomain = 0 ";
      exit(0); 
   }

   // get pointers to the nodes
   
   int Nd1 = connectedExternalNodes(0);  
   int Nd2 = connectedExternalNodes(1);
   
   theNodes[0] = theDomain->getNode(Nd1);
   theNodes[1] = theDomain->getNode(Nd2);  

   if (theNodes[0] == 0)
   {
      opserr << "NLBeamColumn2d::setDomain: Nd1: ";
      opserr << Nd1 << "does not exist in model\n";
      exit(0);
   }

   if (theNodes[1] == 0) 
   {
      opserr << "NLBeamColumn2d::setDomain: Nd2: ";
      opserr << Nd2 << "does not exist in model\n";
      exit(0);
   }

   // call the DomainComponent class method 
   this->DomainComponent::setDomain(theDomain);
    
   // ensure connected nodes have correct number of dof's
   int dofNode1 = theNodes[0]->getNumberDOF();
   int dofNode2 = theNodes[1]->getNumberDOF();
   
   if ((dofNode1 !=3 ) || (dofNode2 != 3))
   {
      opserr << "NLBeamColumn2d::setDomain(): Nd2 or Nd1 incorrect dof ";
      exit(0);
   }
   
   // initialize the transformation
   if (crdTransf->initialize(theNodes[0], theNodes[1]))
   {
      opserr << "NLBeamColumn2d::setDomain(): Error initializing coordinate transformation";  
      exit(0);
   }
    
   // get element length
   double L = crdTransf->getInitialLength();
   if (L == 0.0)
   {
      opserr << "NLBeamColumn2d::setDomain(): Zero element length:" << this->getTag();  
      exit(0);
   }

   if (initialFlag == 0) 
     this->initializeSectionHistoryVariables();
}



int
NLBeamColumn2d::commitState()
{
   int err = 0;
   int i = 0;

   // call element commitState to do any base class stuff
   if ((err = this->Element::commitState()) != 0) {
     opserr << "NLBeamColumn2d::commitState () - failed in base class";
     return err;
   }    

   do {
      vscommit[i] = vs[i];
      err = sections[i++]->commitState();
  
   } while (err == 0 && i < nSections);
   
   if (err)
      return err;
   
   // commit the transformation between coord. systems
   if ((err = crdTransf->commitState()) != 0)
      return err;
      
   // commit the element variables state
   kvcommit = kv;
   Secommit = Se;

   /***
   if (this->getTag() == 75) {
    for (int i=0; i<nSections; i++) {
      opserr << "stress: " << sections[i]->getStressResultant();
      opserr << "strain: " << sections[i]->getSectionDeformation();
      opserr << "flexibility: " << sections[i]->getSectionFlexibility();
    }
   }
   ***/

   //   initialFlag = 0;  fmk - commented out, see what happens to Example3.1.tcl if uncommented
   //                         - i have not a clue why, ask remo if he ever gets in contact with us again!

   return err;
}


int NLBeamColumn2d::revertToLastCommit()
{
   int err;
   int i = 0;
   
   do {
      vs[i] = vscommit[i];
      err = sections[i]->revertToLastCommit();

      sections[i]->setTrialSectionDeformation(vs[i]);      
      Ssr[i] = sections[i]->getStressResultant();
      fs[i]  = sections[i]->getSectionFlexibility();
      
      i++;
   } while (err == 0 && i < nSections);
   
       
   if (err)
      return err;
   
   // revert the transformation to last commit
   if ((err = crdTransf->revertToLastCommit()) != 0)
      return err;
     
   // revert the element state to last commit
   Se   = Secommit;
   kv   = kvcommit;
   
   initialFlag = 0;
   // this->update();

   return err;
}


int NLBeamColumn2d::revertToStart()
{
   // revert the sections state to start
   int err;
   int i = 0;
     
   do {
       fs[i].Zero();
       vs[i].Zero();
       Ssr[i].Zero();
       err = sections[i++]->revertToStart();
 
   }while (err == 0 && i < nSections);

   if (err)
      return err;
   
   // revert the transformation to start
   if ((err = crdTransf->revertToStart()) != 0)
      return err;
  
   // revert the element state to start
   Secommit.Zero();
   kvcommit.Zero();

   Se.Zero();
   kv.Zero();

   initialFlag = 0;
   // this->update();
   return err;
}




const Matrix &
NLBeamColumn2d::getTangentStiff(void)
{
  crdTransf->update();	// Will remove once we clean up the corotational 2d transformation -- MHS
  return crdTransf->getGlobalStiffMatrix(kv, Se);
}


const Matrix &
NLBeamColumn2d::getInitialStiff(void)
{
  
  // check for quick return
  if (Ki != 0)
    return *Ki;

  const Matrix &xi_pt  = quadRule.getIntegrPointCoords(nSections);
  const Vector &weight = quadRule.getIntegrPointWeights(nSections);
  
  static Matrix f(NEBD,NEBD);   // element flexibility matrix
  
  static Matrix I(NEBD,NEBD);   // an identity matrix for matrix inverse
  int i;
  
  I.Zero();
  for (i=0; i<NEBD; i++)
    I(i,i) = 1.0;
  
  double L = crdTransf->getInitialLength();
  double oneOverL  = 1.0/L;
  
  // initialize f and vr for integration
  f.Zero();
  for (i=0; i<nSections; i++) {
	  
    int order      = sections[i]->getOrder();
    const ID &code = sections[i]->getType();
    
    Matrix fb(workArea, order, NEBD);
    
    double xL  = xi_pt(i,0);
    double xL1 = xL-1.0;
    
    // get section flexibility matrix
    const Matrix &fSec = sections[i]->getInitialFlexibility();
	      
    // integrate element flexibility matrix
    // f = f + (b^ fs * b) * weight(i);
    //f.addMatrixTripleProduct(1.0, b[i], fs[i], weight(i));
    int jj, ii;
    fb.Zero();
    double tmp;
    for (ii = 0; ii < order; ii++) {
      switch(code(ii)) {
      case SECTION_RESPONSE_P:
	for (jj = 0; jj < order; jj++)
	  fb(jj,0) += fSec(jj,ii)*weight(i);
	break;
      case SECTION_RESPONSE_MZ:
	for (jj = 0; jj < order; jj++) {
	  tmp = fSec(jj,ii)*weight(i);
	  fb(jj,1) += xL1*tmp;
	  fb(jj,2) += xL*tmp;
	}
	break;
      case SECTION_RESPONSE_VY:
	for (jj = 0; jj < order; jj++) {
	  tmp = oneOverL*fSec(jj,ii)*weight(i);
	  fb(jj,1) += tmp;
	  fb(jj,2) += tmp;
	}
	break;
      default:
	break;
      }
    }
    for (ii = 0; ii < order; ii++) {
      switch (code(ii)) {
      case SECTION_RESPONSE_P:
	for (jj = 0; jj < 3; jj++)
	  f(0,jj) += fb(ii,jj);
	      break;
      case SECTION_RESPONSE_MZ:
	for (jj = 0; jj < 3; jj++) {
	  tmp = fb(ii,jj);
	  f(1,jj) += xL1*tmp;
	  f(2,jj) += xL*tmp;
	}
	break;
      case SECTION_RESPONSE_VY:
	for (jj = 0; jj < 3; jj++) {
	  tmp = oneOverL*fb(ii,jj);
	  f(1,jj) += tmp;
	  f(2,jj) += tmp;
	}
	break;
      default:
	break;
      }
    }
  }  
  
  f  *= L;

  // calculate element stiffness matrix
  // invert3by3Matrix(f, kv);
  static Matrix kvInit(NEBD, NEBD);
  if (f.Solve(I, kvInit) < 0)
    opserr << "NLBeamColumn2d::getInitialStiff() -- could not invert flexibility\n";
			    

  Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(kvInit));
  
  return *Ki;
    
}
    

const Vector &
NLBeamColumn2d::getResistingForce(void)
{
  // Will remove once we clean up the corotational 2d transformation -- MHS
  crdTransf->update();

  Vector p0Vec(p0, 3);

  return crdTransf->getGlobalResistingForce(Se, p0Vec);
}



void
NLBeamColumn2d::initializeSectionHistoryVariables (void)
{
    for (int i = 0; i < nSections; i++)
    {
	int order = sections[i]->getOrder();
	
	fs[i] = Matrix(order,order);
	vs[i] = Vector(order);
        Ssr[i] = Vector(order);

	vscommit[i] = Vector(order);
    }
}


/********* NEWTON , SUBDIVIDE AND INITIAL ITERATIONS ********************
 */
int NLBeamColumn2d::update()
{
  // if have completed a recvSelf() - do a revertToLastCommit
  // to get Ssr, etc. set correctly
  if (initialFlag == 2)
    this->revertToLastCommit();

  // update the transformation
  crdTransf->update();
       
  // get basic displacements and increments
  const Vector &v = crdTransf->getBasicTrialDisp();    
  static Vector dv(NEBD);
  dv = crdTransf->getBasicIncrDeltaDisp();    

  static Vector vin(NEBD);
  vin = v;
  vin -= dv;

  if (initialFlag != 0 && dv.Norm() <= DBL_EPSILON && sp == 0)
    return 0;

  const Matrix &xi_pt  = quadRule.getIntegrPointCoords(nSections);
  const Vector &weight = quadRule.getIntegrPointWeights(nSections);
  
  static Vector vr(NEBD);       // element residual displacements
  static Matrix f(NEBD,NEBD);   // element flexibility matrix
  
  static Matrix I(NEBD,NEBD);   // an identity matrix for matrix inverse
  double dW;                    // section strain energy (work) norm 
  int i, j;
  
  I.Zero();
  for (i=0; i<NEBD; i++)
    I(i,i) = 1.0;

  int numSubdivide = 1;
  bool converged = false;
  static Vector dSe(NEBD);
  static Vector dvToDo(NEBD);
  static Vector dvTrial(NEBD);
  static Vector SeTrial(NEBD);
  static Matrix kvTrial(NEBD, NEBD);

  dvToDo = dv;
  dvTrial = dvToDo;

  static double factor = 10;
  double L = crdTransf->getInitialLength();
  double oneOverL  = 1.0/L;  

  maxSubdivisions = 10;


  // fmk - modification to get compatable ele forces and deformations 
  //   for a change in deformation dV we try first a newton iteration, if
  //   that fails we try an initial flexibility iteration on first iteration 
  //   and then regular newton, if that fails we use the initial flexiblity
  //   for all iterations.
  //
  //   if they both fail we subdivide dV & try to get compatable forces
  //   and deformations. if they work and we have subdivided we apply
  //   the remaining dV.


  while (converged == false && numSubdivide <= maxSubdivisions) {

    // try regular newton (if l==0), or
    // initial tangent on first iteration then regular newton (if l==1), or 
    // initial tangent iterations (if l==2)

    for (int l=0; l<3; l++) {

      //      if (l == 1) l = 2;
      SeTrial = Se;
      kvTrial = kv;
      for (i=0; i<nSections; i++) {
	vsSubdivide[i] = vs[i];
	fsSubdivide[i] = fs[i];
	SsrSubdivide[i] = Ssr[i];
      }

      // calculate nodal force increments and update nodal forces      
      // dSe = kv * dv;
      dSe.addMatrixVector(0.0, kvTrial, dvTrial, 1.0);
      SeTrial += dSe;

      if (initialFlag != 2) {

	int numIters = maxIters;
	if (l == 1) 
	  numIters = 10*maxIters; // allow 10 times more iterations for initial tangent
	
	for (j=0; j <numIters; j++) {
	  
	  // initialize f and vr for integration
	  f.Zero();
	  vr.Zero();
	  
	  for (i=0; i<nSections; i++) {
	    
	    int order      = sections[i]->getOrder();
	    const ID &code = sections[i]->getType();
	    
	    static Vector Ss;
	    static Vector dSs;
	    static Vector dvs;
	    static Matrix fb;
	    
	    Ss.setData(workArea, order);
	    dSs.setData(&workArea[order], order);
	    dvs.setData(&workArea[2*order], order);
	    fb.setData(&workArea[3*order], order, NEBD);
	    
	    double xL  = xi_pt(i,0);
	    double xL1 = xL-1.0;
	    
	    // calculate total section forces
	    // Ss = b*Se + bp*currDistrLoad;
	    // Ss.addMatrixVector(0.0, b[i], Se, 1.0);
	    int ii;
	    for (ii = 0; ii < order; ii++) {
	      switch(code(ii)) {
	      case SECTION_RESPONSE_P:
		Ss(ii) = SeTrial(0);
		break;
	      case SECTION_RESPONSE_MZ:
		Ss(ii) =  xL1*SeTrial(1) + xL*SeTrial(2);
		break;
	      case SECTION_RESPONSE_VY:
		Ss(ii) = oneOverL*(SeTrial(1)+SeTrial(2));
		break;
	      default:
		Ss(ii) = 0.0;
		break;
	      }
	    }
	    
	    // Add the effects of element loads, if present
	    if (sp != 0) {
	      const Matrix &s_p = *sp;
	      for (ii = 0; ii < order; ii++) {
		switch(code(ii)) {
		case SECTION_RESPONSE_P:
		  Ss(ii) += s_p(0,i);
		  break;
		case SECTION_RESPONSE_MZ:
		  Ss(ii) += s_p(1,i);
		  break;
		case SECTION_RESPONSE_VY:
		  Ss(ii) += s_p(2,i);
		  break;
		default:
		  break;
		}
	      }
	    }
	    
	    // dSs = Ss - Ssr[i];
	    dSs = Ss;
	    dSs.addVector(1.0, SsrSubdivide[i], -1.0);
	    

	    // compute section deformation increments
	    if (l == 0) {

	      //  regular newton 
	      //    vs += fs * dSs;     

	      dvs.addMatrixVector(0.0, fsSubdivide[i], dSs, 1.0);
	    } else if (l == 2) {

	      //  newton with initial tangent if first iteration
	      //    vs += fs0 * dSs;     
	      //  otherwise regular newton 
	      //    vs += fs * dSs;     

	      if (j == 0) {
		const Matrix &fs0 = sections[i]->getInitialFlexibility();

		dvs.addMatrixVector(0.0, fs0, dSs, 1.0);
	      } else
		dvs.addMatrixVector(0.0, fsSubdivide[i], dSs, 1.0);

	    } else {
	     
	      //  newton with initial tangent
	      //    vs += fs0 * dSs;     
	    
	      const Matrix &fs0 = sections[i]->getInitialFlexibility();
	      dvs.addMatrixVector(0.0, fs0, dSs, 1.0);
	    }
	    
	    // set section deformations
	    if (initialFlag != 0)
	      vsSubdivide[i] += dvs;
	    
	    sections[i]->setTrialSectionDeformation(vsSubdivide[i]);
	    
	    // get section resisting forces
	    SsrSubdivide[i] = sections[i]->getStressResultant();
	    
	    // get section flexibility matrix
	    fsSubdivide[i] = sections[i]->getSectionFlexibility();
	    
	    // calculate section residual deformations
	    // dvs = fs * (Ss - Ssr);
	    dSs = Ss;
	    dSs.addVector(1.0, SsrSubdivide[i], -1.0);  // dSs = Ss - Ssr[i];
	    
	    dvs.addMatrixVector(0.0, fsSubdivide[i], dSs, 1.0);
	    
	    // integrate element flexibility matrix
	    // f = f + (b^ fs * b) * weight(i);
	    //f.addMatrixTripleProduct(1.0, b[i], fs[i], weight(i));
	    int jj;
	    const Matrix &fSec = fsSubdivide[i];
	    fb.Zero();
	    double tmp;
	    for (ii = 0; ii < order; ii++) {
	      switch(code(ii)) {
	      case SECTION_RESPONSE_P:
		for (jj = 0; jj < order; jj++)
		  fb(jj,0) += fSec(jj,ii)*weight(i);
		break;
	      case SECTION_RESPONSE_MZ:
		for (jj = 0; jj < order; jj++) {
		  tmp = fSec(jj,ii)*weight(i);
		  fb(jj,1) += xL1*tmp;
		  fb(jj,2) += xL*tmp;
		}
		break;
	      case SECTION_RESPONSE_VY:
		for (jj = 0; jj < order; jj++) {
		  tmp = oneOverL*fSec(jj,ii)*weight(i);
		  fb(jj,1) += tmp;
		  fb(jj,2) += tmp;
		}
		break;
	      default:
		break;
	      }
	    }
	    
	    for (ii = 0; ii < order; ii++) {
	      switch (code(ii)) {
	      case SECTION_RESPONSE_P:
		for (jj = 0; jj < 3; jj++)
		  f(0,jj) += fb(ii,jj);
		break;
	      case SECTION_RESPONSE_MZ:
		for (jj = 0; jj < 3; jj++) {
		  tmp = fb(ii,jj);
		  f(1,jj) += xL1*tmp;
		  f(2,jj) += xL*tmp;
		}
		break;
	      case SECTION_RESPONSE_VY:
		for (jj = 0; jj < 3; jj++) {
		  tmp = oneOverL*fb(ii,jj);
		  f(1,jj) += tmp;
		  f(2,jj) += tmp;
		}
		break;
	      default:
		break;
	      }
	    }
	    
	    // integrate residual deformations
	    // vr += (b^ (vs + dvs)) * weight(i);
	    //vr.addMatrixTransposeVector(1.0, b[i], vs[i] + dvs, weight(i));
	    dvs.addVector(1.0, vsSubdivide[i], 1.0);
	    double dei;
	    for (ii = 0; ii < order; ii++) {
	      dei = dvs(ii)*weight(i);
	      switch(code(ii)) {
	      case SECTION_RESPONSE_P:
		vr(0) += dei; break;
	      case SECTION_RESPONSE_MZ:
		vr(1) += xL1*dei; vr(2) += xL*dei; break;
	      case SECTION_RESPONSE_VY:
		tmp = oneOverL*dei;
		vr(1) += tmp; vr(2) += tmp; break;
	      default:
		break;
	      }
	    }
	  }
	  
	  f  *= L;
	  vr *= L;
	  
	  // calculate element stiffness matrix
	  // invert3by3Matrix(f, kv);
	  
	  if (f.Solve(I, kvTrial) < 0)
	    opserr << "NLBeamColumn2d::update() -- could not invert flexibility\n";
				    

	  // dv = vin + dvTrial  - vr
	  dv = vin;
	  dv += dvTrial;
	  dv -= vr;
	  
	  // dv.addVector(1.0, vr, -1.0);
	  
	  // dSe = kv * dv;
	  dSe.addMatrixVector(0.0, kvTrial, dv, 1.0);

	  dW = dv ^ dSe; 
	  
	  SeTrial += dSe;
	  
	  // check for convergence of this interval
	  if (fabs(dW) < tol) { 
	    
	    // set the target displacement
	    dvToDo -= dvTrial;
	    vin += dvTrial;
	    
	    // check if we have got to where we wanted
	    if (dvToDo.Norm() <= DBL_EPSILON) {
	      converged = true;
	      
	    } else {  // we convreged but we have more to do

	      // reset variables for start of next subdivision
	      dvTrial = dvToDo;
	      numSubdivide = 1;  // NOTE setting subdivide to 1 again maybe too much
	    }
	    
	    // set kv, vs and Se values
	    kv = kvTrial;
	    Se = SeTrial;
	    
	    for (int k=0; k<nSections; k++) {
	      vs[k] = vsSubdivide[k];
	      fs[k] = fsSubdivide[k];
	      Ssr[k] = SsrSubdivide[k];
	    }

	    // break out of j & l loops
	    j = numIters+1;
	    l = 4;

	  } else {   //  if (fabs(dW) < tol) { 

	    // if we have failed to convrege for all of our newton schemes
	    // - reduce step size by the factor specified

	    if (j == (numIters-1) && (l == 2)) {
	      dvTrial /= factor;
	      numSubdivide++;
	    }
	  }

	} // for (j=0; j<numIters; j++)
      } // if (initialFlag != 2)
    } // for (int l=0; l<2; l++)
  } // while (converged == false)


  // if fail to converge we return an error flag & print an error message

  if (converged == false) {
    opserr << "WARNING - NLBeamColumn2d::update - failed to get compatable ";
    opserr << "element forces & deformations for element: ";
    opserr << this->getTag() << "(dW: << " << dW << ")\n";
    return -1;
  }

  initialFlag = 1;

  return 0;
}


void 
NLBeamColumn2d::getGlobalDispls(Vector &dg) const
{
   // determine global displacements
   const Vector &disp1 = theNodes[0]->getTrialDisp();
   const Vector &disp2 = theNodes[1]->getTrialDisp();

   for (int i = 0; i < NND; i++)
   {
      dg(i)     = disp1(i);
      dg(i+NND) = disp2(i);
   }
}



void NLBeamColumn2d::getGlobalAccels(Vector &ag) const
{
   // determine global displacements
   const Vector &accel1 = theNodes[0]->getTrialAccel();
   const Vector &accel2 = theNodes[1]->getTrialAccel();

   for (int i = 0; i < NND; i++)
   {
      ag(i)     = accel1(i);
      ag(i+NND) = accel2(i);
   }
}



void NLBeamColumn2d::getForceInterpolatMatrix(double xi, Matrix &b, const ID &code)
{
   b.Zero();

   double L = crdTransf->getInitialLength();
   for (int i = 0; i < code.Size(); i++)
   {
      switch (code(i))
      {
         case SECTION_RESPONSE_MZ:		// Moment, Mz, interpolation
	    b(i,1) = xi - 1.0;
	    b(i,2) = xi;
	    break;
	 case SECTION_RESPONSE_P:		// Axial, P, interpolation
	    b(i,0) = 1.0;
	    break;
	 case SECTION_RESPONSE_VY:		// Shear, Vy, interpolation
	    b(i,1) = b(i,2) = 1.0/L;
	    break;
	 default:
	    break;
      }
   }
}


void NLBeamColumn2d::getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code)
{
   bp.Zero();

   double L = crdTransf->getInitialLength();
   for (int i = 0; i < code.Size(); i++)
   {
      switch (code(i))
      {
         case SECTION_RESPONSE_MZ:		// Moment, Mz, interpolation
	    bp(i,1) = xi*(xi-1)*L*L/2;
	    break;
	 case SECTION_RESPONSE_P:		// Axial, P, interpolation
	    bp(i,0) = (1-xi)*L;
	    break;
	 case SECTION_RESPONSE_VY:		// Shear, Vy, interpolation
	    bp(i,1) = (xi-0.5)*L;
	    break;
	 default:
	    break;
      }
   }
}


const Matrix &
NLBeamColumn2d::getMass(void)
{ 
  theMatrix.Zero();

  double L = crdTransf->getInitialLength();
  if (rho != 0.0)
    theMatrix(0,0) = theMatrix(1,1) = theMatrix(3,3) = theMatrix(4,4) = 0.5*L*rho;
   
    return theMatrix;
}



void 
NLBeamColumn2d::zeroLoad(void)
{
  if (sp != 0) {
    sp->Zero();

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
  }

  load.Zero();
}

int
NLBeamColumn2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  
  if (sp == 0) {
    sp = new Matrix(3,nSections);
    if (sp == 0) {
      opserr << "NLBeamColumn2d::addLoad -- out of memory\n";
      exit(-1);
    }
  }

  const Matrix &xi_pt = quadRule.getIntegrPointCoords(nSections);
  double L = crdTransf->getInitialLength();

  if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wa = data(1)*loadFactor;  // Axial
    double wy = data(0)*loadFactor;  // Transverse

    Matrix &s_p = *sp;

    // Accumulate applied section forces due to element loads
    for (int i = 0; i < nSections; i++) {
      double x = xi_pt(i,0)*L;
      // Axial
      s_p(0,i) += wa*(L-x);
      // Moment
      s_p(1,i) += wy*0.5*x*(x-L);
      // Shear
      s_p(2,i) += wy*(x-0.5*L);
    }

    // Accumulate reactions in basic system
    p0[0] -= wa*L;
    double V = 0.5*wy*L;
    p0[1] -= V;
    p0[2] -= V;
  }
  else if (type == LOAD_TAG_Beam2dPointLoad) {
    double P = data(0)*loadFactor;
    double N = data(1)*loadFactor;
    double aOverL = data(2);
    double a = aOverL*L;

    double V1 = P*(1.0-aOverL);
    double V2 = P*aOverL;

    Matrix &s_p = *sp;

    // Accumulate applied section forces due to element loads
    for (int i = 0; i < nSections; i++) {
      double x = xi_pt(i,0)*L;
      if (x <= a) {
	s_p(0,i) += N;
	s_p(1,i) -= x*V1;
	s_p(2,i) -= V1;
      }
      else {
	s_p(1,i) -= (L-x)*V2;
	s_p(2,i) += V2;
      }
    }

    // Accumulate reactions in basic system
    p0[0] -= N;
    p0[1] -= V1;
    p0[2] -= V2;
  }

  else {
    opserr << "NLBeamColumn2d::addLoad() -- load type unknown for element with tag: " << this->getTag() << endln;
    return -1;
  }
  
  return 0;
}


int 
NLBeamColumn2d::addInertiaLoadToUnbalance(const Vector &accel)
{
  // Check for a quick return
  if (rho == 0.0)
    return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);    

  double L = crdTransf->getInitialLength();
  double m = 0.5*rho*L;

  load(0) -= m*Raccel1(0);
  load(1) -= m*Raccel1(1);
  load(3) -= m*Raccel2(0);
  load(4) -= m*Raccel2(1);

  return 0;
}

const Vector &
NLBeamColumn2d::getResistingForceIncInertia()
{	
  // Compute the current resisting force
  theVector = this->getResistingForce();

  // Check for a quick return
  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
    
    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;
    
    theVector(0) += m*accel1(0);
    theVector(1) += m*accel1(1);
    theVector(3) += m*accel2(0);
    theVector(4) += m*accel2(1);
    
    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      theVector += this->getRayleighDampingForces();

  } else {
    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      theVector += this->getRayleighDampingForces();
  }

  return theVector;
}



bool
NLBeamColumn2d::isSubdomain(void)
{
    return false;
}



int
NLBeamColumn2d::sendSelf(int commitTag, Channel &theChannel)
{  
  // place the integer data into an ID
  int dbTag = this->getDbTag();
  int i, j , k;
  int loc = 0;
  
  static ID idData(9);  // one bigger than needed so no clash later
  idData(0) = this->getTag();
  idData(1) = connectedExternalNodes(0);
  idData(2) = connectedExternalNodes(1);
  idData(3) = nSections;
  idData(4) = maxIters;
  idData(5) = initialFlag;
  idData(6) = crdTransf->getClassTag();
  int crdTransfDbTag  = crdTransf->getDbTag();
  if (crdTransfDbTag  == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag  != 0) 
      crdTransf->setDbTag(crdTransfDbTag);
  }
  idData(7) = crdTransfDbTag;
  

  if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
    opserr << "NLBeamColumn2d::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  // send the coordinate transformation
  
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
    opserr << "NLBeamColumn2d::sendSelf() - failed to send crdTranf\n";
    return -1;
  }      

  
  //
  // send an ID for the sections containing each sections dbTag and classTag
  // if section ha no dbTag get one and assign it
  //

  ID idSections(2*nSections);
  loc = 0;
  for (i = 0; i<nSections; i++) {
    int sectClassTag = sections[i]->getClassTag();
    int sectDbTag = sections[i]->getDbTag();
    if (sectDbTag == 0) {
      sectDbTag = theChannel.getDbTag();
      sections[i]->setDbTag(sectDbTag);
    }

    idSections(loc) = sectClassTag;
    idSections(loc+1) = sectDbTag;
    loc += 2;
  }

  if (theChannel.sendID(dbTag, commitTag, idSections) < 0)  {
    opserr << "NLBeamColumn2d::sendSelf() - failed to send ID data\n";
			    
    return -1;
  }    

  //
  // send the sections
  //
  
  for (j = 0; j<nSections; j++) {
    if (sections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "NLBeamColumn2d::sendSelf() - section: " << j << " failed to send itself\n";
      return -1;
    }
  }
  
  // into a vector place distrLoadCommit, rho, UeCommit, Secommit and kvcommit
  int secDefSize = 0;
  for (i = 0; i < nSections; i++) {
     int size = sections[i]->getOrder();
     secDefSize   += size;
  }

  Vector dData(1+1+NEBD+NEBD*NEBD+secDefSize); 
  loc = 0;

  // place double variables into Vector
  dData(loc++) = rho;
  dData(loc++) = tol;
  
  // put  distrLoadCommit into the Vector
  //  for (i=0; i<NL; i++) 
  //dData(loc++) = distrLoadcommit(i);

  // place kvcommit into vector
  for (i=0; i<NEBD; i++) 
    dData(loc++) = Secommit(i);

  // place kvcommit into vector
  for (i=0; i<NEBD; i++) 
     for (j=0; j<NEBD; j++)
        dData(loc++) = kvcommit(i,j);
  
  // place vscommit into vector
  for (k=0; k<nSections; k++)
     for (i=0; i<sections[k]->getOrder(); i++)
	dData(loc++) = (vscommit[k])(i);
  
  if (theChannel.sendVector(dbTag, commitTag, dData) < 0) {
     opserr << "NLBeamColumn2d::sendSelf() - failed to send Vector data\n";
	 		     
     return -1;
  }    

  return 0;
}    


int
NLBeamColumn2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i,j,k;
  
  static ID idData(9); // one bigger than needed 

  if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
    opserr << "NLBeamColumn2d::recvSelf() - failed to recv ID data\n";
    return -1;
  }    

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);
  maxIters = idData(4);
  initialFlag = idData(5);
  
  int crdTransfClassTag = idData(6);
  int crdTransfDbTag = idData(7);

  // create a new crdTransf object if one needed
  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
      if (crdTransf != 0)
	  delete crdTransf;

      crdTransf = theBroker.getNewCrdTransf2d(crdTransfClassTag);

      if (crdTransf == 0) {
	opserr << "NLBeamColumn2d::recvSelf() - failed to obtain a CrdTrans object with classTag " << 
	  crdTransfClassTag << endln;
	return -2;	  
      }
  }

  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf obkject
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0)  
  {
     opserr << "NLBeamColumn2d::sendSelf() - failed to recv crdTranf\n";
     return -3;
  }      
  
  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*idData(3));
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "NLBeamColumn2d::recvSelf() - failed to recv ID data\n";
    return -1;
  }    

  //
  // now receive the sections
  //
  
  if (nSections != idData(3)) {

    //
    // we do not have correct number of sections, must delete the old and create
    // new ones before can recvSelf on the sections
    //

    // delete the old
    if (nSections != 0) {
      for (int i=0; i<nSections; i++)
	delete sections[i];
      delete [] sections;
    }

    // create a section and recvSelf on it
    nSections = idData(3);

    // Delete the old
    if (vscommit != 0)
      delete [] vscommit;

    // Allocate the right number
    vscommit = new Vector[nSections];
    if (vscommit == 0) {
      opserr << "NLBeamColumn2d::recvSelf -- failed to allocate vscommit array\n";
      return -1;
    }

    // Delete the old
    if (fs != 0)
      delete [] fs;

    // Allocate the right number
    fs = new Matrix[nSections];  
    if (fs == 0) {
      opserr << "NLBeamColumn2d::recvSelf -- failed to allocate fs array\n";
      return -1;
    }

    // Delete the old
    if (vs != 0)
      delete [] vs;

    // Allocate the right number
    vs = new Vector[nSections];  
    if (vs == 0) {
      opserr << "NLBeamColumn2d::recvSelf -- failed to allocate vs array\n";
      return -1;
    }

    // Delete the old
    if (Ssr != 0)
      delete [] Ssr;
    
    // Allocate the right number
    Ssr = new Vector[nSections];  
    if (Ssr == 0) {
      opserr << "NLBeamColumn2d::recvSelf -- failed to allocate Ssr array\n";
      return -1;
    }

    // create a new array to hold pointers
    sections = new SectionForceDeformation *[idData(3)];
    if (sections == 0) {
      opserr << "NLBeamColumn2d::recvSelf() - " <<
	"out of memory creating sections array of size" << idData(3) << endln;
      exit(-1);
    }    

    loc = 0;
    
    for (i=0; i<nSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;
      sections[i] = theBroker.getNewSection(sectClassTag);
      if (sections[i] == 0) {
	opserr << "NLBeamColumn2d::recvSelf() - " <<
	  "Broker could not create Section of class type" << sectClassTag << endln;
	exit(-1);
      }
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "NLBeamColumn2d::recvSelf() - section " << i << "failed to recv itself\n";
	return -1;
      }     
    }

    this->initializeSectionHistoryVariables();

  } else {

    // 
    // for each existing section, check it is of correct type
    // (if not delete old & create a new one) then recvSelf on it
    //
    
    loc = 0;
    for (i=0; i<nSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;

      // check of correct type
      if (sections[i]->getClassTag() !=  sectClassTag) {
	// delete the old section[i] and create a new one
	delete sections[i];
	sections[i] = theBroker.getNewSection(sectClassTag);
	if (sections[i] == 0) {
	  opserr << "NLBeamColumn2d::recvSelf() - " << 
	    "Broker could not create Section of class type" <<sectClassTag << endln;
	  exit(-1);
	}
      }

      // recvvSelf on it
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "NLBeamColumn2d::recvSelf() - section " << 
	  i << "failed to recv itself\n";
	return -1;
      }     
    }
  }
  
  // into a vector place distrLoadCommit, rho, UeCommit, Secommit and kvcommit
  int secDefSize = 0;
  for (int ii = 0; ii < nSections; ii++) {
     int size = sections[ii]->getOrder();
     secDefSize   += size;
  }
  
  Vector dData(1+1+NEBD+NEBD*NEBD+secDefSize);   
  
  if (theChannel.recvVector(dbTag, commitTag, dData) < 0)  {
    opserr << "NLBeamColumn2d::sendSelf() - failed to send Vector data\n";
			    
    return -1;
  }    
  
  loc = 0;

  // place double variables into Vector
  rho = dData(loc++);
  tol = dData(loc++);

  // put  distrLoadCommit into the Vector
  //for (i=0; i<NL; i++) 
  // distrLoad(i) = dData(loc++);

  // place kvcommit into vector
  for (i=0; i<NEBD; i++) 
    Secommit(i) = dData(loc++);

  // place kvcommit into vector
  for (i=0; i<NEBD; i++) 
     for (j=0; j<NEBD; j++)
        kvcommit(i,j) = dData(loc++);

  kv   = kvcommit;
  Se   = Secommit;

  for (k = 0; k < nSections; k++) {
    int order = sections[k]->getOrder();

    // place vscommit into vector
    vscommit[k] = Vector(order);
    for (i = 0; i < order; i++)
      (vscommit[k])(i) = dData(loc++);
  }

  initialFlag = 2;  

  return 0;
}




void NLBeamColumn2d::compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const
{

   // update the transformation
   crdTransf->update();
       
   // get basic displacements and increments
   static Vector ub(NEBD);
   ub = crdTransf->getBasicTrialDisp();    
  
   // get integration point positions and weights
   const Matrix &xi_pt  = quadRule.getIntegrPointCoords(nSections);

   // setup Vandermode and CBDI influence matrices
   int i;
   double xi;
 
   // get CBDI influence matrix
   Matrix ls(nSections, nSections);
   double L = crdTransf->getInitialLength();
   getCBDIinfluenceMatrix(nSections, xi_pt, L, ls);

   // get section curvatures
   Vector kappa(nSections);  // curvature
   Vector vs;              // section deformations 

   for (i=0; i<nSections; i++)
   {
       // THIS IS VERY INEFFICIENT ... CAN CHANGE LATER
       int sectionKey = 0;
       const ID &code = sections[i]->getType();
       int ii;
       for (ii = 0; ii < code.Size(); ii++)
	   if (code(ii) == SECTION_RESPONSE_MZ)
	   {
	       sectionKey = ii;
	       break;
	   }

       if (ii == code.Size()) {
	 opserr << "FATAL NLBeamColumn2d::compSectionDispls - section does not provide Mz response\n";
	 exit(-1);
       }
			
       // get section deformations
       vs = sections[i]->getSectionDeformation();
       kappa(i) = vs(sectionKey);
   }

       //opserr << "kappa: " << kappa;   

   Vector w(nSections);
   static Vector xl(NDM), uxb(NDM);
   static Vector xg(NDM), uxg(NDM); 

   // w = ls * kappa;  
   w.addMatrixVector (0.0, ls, kappa, 1.0);
   
   for (i=0; i<nSections; i++)
   {
      xi = xi_pt(i,0);

      xl(0) = xi * L;
      xl(1) = 0;

      // get section global coordinates
      sectionCoords[i] = crdTransf->getPointGlobalCoordFromLocal(xl);

      // compute section displacements
      uxb(0) = xi * ub(0); // consider linear variation for axial displacement. CHANGE LATER!!!!!!!!!!
      uxb(1) = w(i);
             
      // get section displacements in global system 
      sectionDispls[i] = crdTransf->getPointGlobalDisplFromBasic(xi, uxb);
   }	       
}

   



void
NLBeamColumn2d::Print(OPS_Stream &s, int flag)
{
   if (flag == -1) { 
    int eleTag = this->getTag();
    s << "NL_BEAM\t" << eleTag << "\t";
    s << sections[0]->getTag() << "\t" << sections[nSections-1]->getTag();
    s  << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
    s << "0\t0.0000000\n";
   } 

   else if (flag == 1) { 
      s << "\nElement: " << this->getTag() << " Type: NLBeamColumn2d ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << nSections;
      s << "\tMass density: " << rho;
      for (int i = 0; i < nSections; i++)
         s << "\nSection "<<i<<" :" << *sections[i];
   } 
   
   else {
      s << "\nElement: " << this->getTag() << " Type: NLBeamColumn2d ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << nSections;
      s << "\tMass density: " << rho << endln;
      double P  = Secommit(0);
      double M1 = Secommit(1);
      double M2 = Secommit(2);
      double L = crdTransf->getInitialLength();
      double V = (Secommit(1)+Secommit(2))/L;
      theVector(1) = V;
      theVector(4) = -V;
      s << "\tEnd 1 Forces (P V M): " << -P+p0[0] << " "
	<< V+p0[1] << " " << M1 << endln;
      s << "\tEnd 2 Forces (P V M): " << P << " "
	<< -V+p0[2] << " " << M2 << endln;
   }
}


OPS_Stream &operator<<(OPS_Stream &s, NLBeamColumn2d &E)
{
    E.Print(s);
    return s;
}



int
NLBeamColumn2d::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
  if (displayMode == 2) {
    // first determine the end points of the beam based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	

    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();

    static Vector v1(3);
    static Vector v2(3);
    
    for (int i = 0; i < 2; i++) {
      v1(i) = end1Crd(i) + end1Disp(i)*fact;
      v2(i) = end2Crd(i) + end2Disp(i)*fact;    
    }
    
    return theViewer.drawLine (v1, v2, 1.0, 1.0);
  }
  
  else if (displayMode == 1) {
    // first determine the two end points of the element based on
    //  the display factor (a measure of the distorted image)
    
    static Vector v1(3), v2(3);
    
    const Vector &node1Crd = theNodes[0]->getCrds();
    const Vector &node2Crd = theNodes[1]->getCrds();	
    const Vector &node1Disp = theNodes[0]->getDisp();
    const Vector &node2Disp = theNodes[1]->getDisp();    
    
    v1(2) = 0.0;
    v2(2) = 0.0;
              
    
    int i;
    
    // allocate array of vectors to store section coordinates and displacements
    static Vector *displs = 0;
    static Vector *coords = 0;
    static int numCoords = 0;
    if (numCoords < maxNumSections) {
      if (coords != 0) {
	delete [] coords;
	delete [] displs;
      }
      coords = new Vector [maxNumSections];
      displs = new Vector [maxNumSections];
      if (coords == 0 || displs == 0) {
	opserr << "NLBeamColumn2d::displaySelf() -- failed to allocate coords array";   
	exit(-1);
      }
      numCoords = maxNumSections;
      for (i = 0; i < nSections; i++)
	coords[i] = Vector(2);
      
      for (i = 0; i < nSections; i++)
	displs[i] = Vector(2);
    }
    
    int error;
    
    this->compSectionDisplacements(coords, displs);
    
    v1(0) = node1Crd(0) + node1Disp(0)*fact;
    v1(1) = node1Crd(1) + node1Disp(1)*fact;
    
    ///opserr << "v1: " << v1;
    
    // get global displacements and coordinates of each section          
    
    for (i=0; i<nSections; i++) {
      
      v2(0) = (coords[i])(0) + ((displs[i])(0))*fact;
      v2(1) = (coords[i])(1) + ((displs[i])(1))*fact;
      
      error = theViewer.drawLine(v1, v2, 1.0, 1.0);
      
      if (error)
	return error;
      v1 = v2;
      
    }
    
    v2(0) = node2Crd(0) + node2Disp(0)*fact;
    v2(1) = node2Crd(1) + node2Disp(1)*fact;
    
    error = theViewer.drawLine(v1, v2, 1.0, 1.0);
    
    if (error)
      return error;
  } 
  
  else if (displayMode == 3) {
    
    
    // plot the curvatures
    // first determine the two end points of the element based on
    //  the display factor (a measure of the distorted image)
    
    static Vector v1(NDM), v2(NDM);
    
    const Vector &node1Crd = theNodes[0]->getCrds();
    const Vector &node2Crd = theNodes[1]->getCrds();	
    
    v1(2) = 0.0;
    v2(2) = 0.0;
    
    // subdivide element into smaller parts to draw the deformed shape
    double x_i, y_i, x_j, y_j, xg_xi0, yg_xi0, xi;
    x_i = node1Crd(0);
    y_i = node1Crd(1);
    x_j = node2Crd(0);
    y_j = node2Crd(1);
    
    // determine displaced position of node i      
    xg_xi0 = x_i;
    yg_xi0 = y_i;
    
    // get integration point positions and weights
    const Matrix &xi_pt = quadRule.getIntegrPointCoords(nSections);
    
    // get section curvatures
    Vector kappa(nSections); // curvature
    Vector vs; // section deformations 
    int i;
    for (i=0; i<nSections; i++) {
      // THIS IS VERY INEFFICIENT ... CAN CHANGE IF RUNS TOO SLOW
      int sectionKey = 0;
      const ID &code = sections[i]->getType();
      int ii;
      for (ii = 0; ii < code.Size(); ii++)
	if (code(ii) == SECTION_RESPONSE_MZ) {
	  sectionKey = ii;
	  break;
	}
      
      if (ii == code.Size()) {
	opserr << "FATAL NLBeamColumn2d::displaySelf - section does not provide Mz response\n";
	exit(-1);
      }
      // get section deformations
      vs = sections[i]->getSectionDeformation();
	 kappa(i) = vs(sectionKey);
       }
       
       double xl_xi, yl_xi, xg_xi, yg_xi;
       int error;
       
       double L = crdTransf->getInitialLength();
       for (i = 0; i< nSections; i++) {
	 xi = xi_pt(i,0);
	 
	 // determine displaced local coordinates of the point xi
	 xl_xi = L * xi;
	 yl_xi = kappa(i) * fact;
	 
	 // rotate to global coordinates
	 xg_xi = cosTheta * xl_xi - sinTheta * yl_xi;
	 yg_xi = sinTheta * xl_xi + cosTheta * yl_xi;
	 
	 // translate to global coordinates
	 xg_xi = xg_xi + x_i;
	 yg_xi = yg_xi + y_i;
	
	 // draw the displaced position of this line segment
	 v1(0) = xg_xi0;
	 v1(1) = yg_xi0;
	 
	 v2(0) = xg_xi;
	 v2(1) = yg_xi;
	 
	 error =  theViewer.drawLine(v1, v2, 1.0, 1.0);	
	 
	 if (error)
	   return error;
	 
	 xg_xi0 = xg_xi;
	 yg_xi0 = yg_xi;
       }

      v1(0) = xg_xi0;
      v1(1) = yg_xi0;
      v2(0) = x_j;
      v2(1) = y_j;

      error =  theViewer.drawLine(v1, v2, 1.0, 1.0);	
	       
      return error;
   
   }
   return 0;
}

Response*
NLBeamColumn2d::setResponse(const char **argv, int argc, Information &eleInformation)
{
    //
    // we compare argv[0] for known response types 
    //

    // global force - 
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
	|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
		return new ElementResponse(this, 1, theVector);

    // local force -
    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
	return new ElementResponse(this, 2, theVector);
    
    // section response -
    else if (strcmp(argv[0],"section") ==0) {
      if (argc <= 2)
	return 0;
	
      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= nSections)
	return sections[sectionNum-1]->setResponse(&argv[2], argc-2, eleInformation);
      else
	return 0;
    }
    
    else
      return 0;
}

int 
NLBeamColumn2d::getResponse(int responseID, Information &eleInfo)
{
  double V;

  switch (responseID) {
    case 1:  // global forces
      return eleInfo.setVector(this->getResistingForce());

    case 2:
      theVector(3) =  Se(0);
      theVector(0) = -Se(0)+p0[0];
      theVector(2) = Se(1);
      theVector(5) = Se(2);
      V = (Se(1)+Se(2))/crdTransf->getInitialLength();
      theVector(1) =  V+p0[1];
      theVector(4) = -V+p0[2];
      return eleInfo.setVector(theVector);
      
    default: 
	  return -1;
  }
}

int
NLBeamColumn2d::setParameter (const char **argv, int argc, Information &info)
{
	//
	// From the parameterID value it should be possible to extract
	// information about:
	//  1) Which parameter is in question. The parameter could
	//     be at element, section, or material level. 
	//  2) Which section and material number (tag) it belongs to. 
	//
	// To accomplish this the parameterID is given the following value:
	//     parameterID = type + 1000*matrTag + 100000*sectionTag
	// ...where 'type' is an integer in the range (1-99) and added 100
	// for each level (from material to section to element). 
	//
	// Example:
	//    If 'E0' (case 2) is random in material #3 of section #5
	//    the value of the parameterID at this (element) level would be:
	//    parameterID = 2 + 1000*3 + 100000*5 = 503002
	//    As seen, all given information can be extracted from this number. 
	//

	// Initial declarations
	int parameterID = 0;

	// If the parameter belongs to the element itself
	if (strcmp(argv[0],"rho") == 0) {
		info.theType = DoubleType;
		return 1;
	}

	// If the parameter is belonging to a section or lower
	else if (strcmp(argv[0],"section") == 0) {

		// For now, no parameters of the section itself:
		if (argc<5) {
			opserr << "For now: cannot handle parameters of the section itself." << endln;
			return -1;
		}

		// Get section and material tag numbers from user input
		int paramSectionTag = atoi(argv[1]);

		// Find the right section and call its setParameter method
		for (int i=0; i<nSections; i++) {
			if (paramSectionTag == sections[i]->getTag()) {
				parameterID = sections[i]->setParameter(&argv[2], argc-2, info);
			}
		}
		
		// Check if the parameterID is valid
		if (parameterID < 0) {
			opserr << "NLBeamColumn2d::setParameter() - could not set parameter. " << endln;
			return -1;
		}
		else {
			// Return the parameterID value (according to the above comments)
			return parameterID;
		}
	}
    
	// Otherwise parameter is unknown for this class
	else {
		return -1;
	}
}

int
NLBeamColumn2d::updateParameter (int parameterID, Information &info)
{
	// If the parameterID value is not equal to 1 it belongs 
	// to section or material further down in the hierarchy. 

	if (parameterID == 1) {

		this->rho = info.theDouble;
		return 0;

	}
	else if (parameterID > 0 ) {

		// Extract the section number
		int sectionNumber = (int)( floor((double)parameterID) / (100000) );

		int ok = -1;
		for (int i=0; i<nSections; i++) {
			if (sectionNumber == sections[i]->getTag()) {
				ok = sections[i]->updateParameter(parameterID, info);
			}
		}

		if (ok < 0) {
			opserr << "NLBeamColumn2d::updateParameter() - could not update parameter. " << endln;
			return ok;
		}
		else {
			return ok;
		}
	}
	else {
		opserr << "NLBeamColumn2d::updateParameter() - could not update parameter. " << endln;
		return -1;
	}       
}
