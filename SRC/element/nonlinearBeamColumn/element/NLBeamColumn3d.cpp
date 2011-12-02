
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
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/element/NLBeamColumn3d.cpp,v $
                                                                        
                                                                        
// File: ~/element/NLBeamColumn3d.C
//
// Written: Remo Magalhaes de Souza (rmsouza.ce.berkeley.edu) on 03/99 
// Revised: rms 06/99 (mass matrix)
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
// Purpose: This file contains the implementation for the NLBeamColumn3d class.
//          NLBeamColumn3d is a materially nonlinear flexibility based frame element.

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iomanip.h>

#include <Information.h>
#include <NLBeamColumn3d.h>
#include <MatrixUtil.h>
#include <GaussQuadRule1d01.h>
#include <GaussLobattoQuadRule1d01.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <G3Globals.h>

#define  NDM   3         // dimension of the problem (3d)
#define  NL    3         // size of uniform load vector
#define  NND   6         // number of nodal dof's
#define  NEGD 12         // number of element global dof's
#define  NEBD  6         // number of element dof's in the basic system

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
NLBeamColumn3d::NLBeamColumn3d():

Element(0,ELE_TAG_NLBeamColumn3d), connectedExternalNodes(2), 
nSections(0), sections(0), node1Ptr(0), node2Ptr(0),
rho(0), maxIters(0), tol(0), initialFlag(0), prevDistrLoad(NL),
K(NEGD,NEGD), m(NEGD,NEGD), d(NEGD,NEGD), P(NEGD), Pinert(NEGD), load(NEGD), Uepr(NEGD),
kv(NEBD,NEBD), Se(NEBD), 
distrLoadcommit(NL), Uecommit(NEGD), kvcommit(NEBD,NEBD), Secommit(NEBD), b(0), bp(0),
fs(0), vs(0), Ssr(0), vscommit(0), crdTransf(0)
{

}

// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points. 
// allocates the necessary space needed by each object
NLBeamColumn3d::NLBeamColumn3d (int tag, int nodeI, int nodeJ, 
                                int numSections, SectionForceDeformation *sectionPtrs[],
                                CrdTransf3d &coordTransf, double massDensPerUnitLength,
                                int maxNumIters, double tolerance):

Element(tag,ELE_TAG_NLBeamColumn3d), connectedExternalNodes(2), 
nSections(numSections), node1Ptr(0), node2Ptr(0),
rho(massDensPerUnitLength), maxIters(maxNumIters), tol(tolerance), 
initialFlag(0), prevDistrLoad(NL), 
K(NEGD,NEGD), m(NEGD,NEGD), d(NEGD,NEGD), P(NEGD), Pinert(NEGD), load(NEGD), 
Uepr(NEGD), kv(NEBD,NEBD), Se(NEBD),  
distrLoadcommit(NL), Uecommit(NEGD), kvcommit(NEBD,NEBD), Secommit(NEBD), b(0), bp(0),
fs(0), vs(0), Ssr(0), vscommit(0)
{
   connectedExternalNodes(0) = nodeI;
   connectedExternalNodes(1) = nodeJ;    

   // get copy of the sections
   
   if (!sectionPtrs)
   {
       cerr << "Error: NLBeamColumn3d::NLBeamColumn3d:  invalid section pointer ";
       exit(-1);
   }	  
   
   sections = new SectionForceDeformation *[nSections];
   if (!sections)
   {
       cerr << "Error: NLBeamColumn3d::NLBeamColumn3d: could not alocate section pointer";
       exit(-1);
   }  
   
   for (int i = 0; i < nSections; i++)
   {
      if (!sectionPtrs[i])
      {
	  cerr << "Error: NLBeamColumn3d::NLBeamColumn3d: section pointer " << i << endl;
          exit(-1);
      }  
       
      sections[i] = sectionPtrs[i]->getCopy();
      if (!sections[i])
      {
	  cerr << "Error: NLBeamColumn3d::NLBeamColumn3d: could not create copy of section " << i << endl;
          exit(-1);
      }
   }

   // get copy of the transformation object   
   crdTransf = coordTransf.getCopy(); 
   if (!crdTransf)
   {
      cerr << "Error: NLBeamColumn3d::NLBeamColumn3d: could not create copy of coordinate transformation object" << endl;
      exit(-1);
   }

   // alocate force interpolation matrices
   
   b  = new Matrix [nSections];
   if (!b)
   {
       cerr << "NLBeamColumn3d::NLBeamColumn3d() -- failed to allocate b array";
       exit(-1);
   }
   
   bp = new Matrix [nSections];
   if (!bp)
   {
       cerr << "NLBeamColumn3d::NLBeamColumn3d() -- failed to allocate bp array";
       exit(-1);
   }

   // alocate section flexibility matrices and section deformation vectors
   fs  = new Matrix [nSections];
   if (!fs)
   {
       cerr << "NLBeamColumn3d::NLBeamColumn3d() -- failed to allocate fs array";
       exit(-1);
   }
   
   vs = new Vector [nSections];
   if (!vs)
   {
       cerr << "NLBeamColumn3d::NLBeamColumn3d() -- failed to allocate vs array";
       exit(-1);
   }

   Ssr = new Vector [nSections];
   if (!Ssr)
   {
       cerr << "NLBeamColumn3d::NLBeamColumn3d() -- failed to allocate Ssr array";
       exit(-1);
   }
 
   vscommit = new Vector [nSections];
   if (!vscommit)
   {
       cerr << "NLBeamColumn3d::NLBeamColumn3d() -- failed to allocate vscommit array";   
       exit(-1);
   }
}


// ~NLBeamColumn3d():
// 	destructor
//      delete must be invoked on any objects created by the object
NLBeamColumn3d::~NLBeamColumn3d()
{
   int i;
   
   if (sections)
   {
      for (i=0; i < nSections; i++)
         if (sections[i])
            delete sections[i];
      delete [] sections;
   }  
   
   if (b)
       delete [] b;
   if (bp)
       delete [] bp;

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
}



int
NLBeamColumn3d::getNumExternalNodes(void) const
{
   return connectedExternalNodes.Size();
}


const ID &
NLBeamColumn3d::getExternalNodes(void)
{
   return connectedExternalNodes;
}


int
NLBeamColumn3d::getNumDOF(void) 
{
   return NEGD;
}

void
NLBeamColumn3d::setDomain(Domain *theDomain)
{
   // check Domain is not null - invoked when object removed from a domain
   if (theDomain == 0)
   {
      node1Ptr = 0;
      node2Ptr = 0;
	  return;
   }

   // get pointers to the nodes
   
   int Nd1 = connectedExternalNodes(0);  
   int Nd2 = connectedExternalNodes(1);
   
   node1Ptr = theDomain->getNode(Nd1);
   node2Ptr = theDomain->getNode(Nd2);  

   if (node1Ptr == 0)
   {
      cerr << "NLBeamColumn3d::setDomain: Nd1: ";
      cerr << Nd1 << "does not exist in model\n";
      exit(0);
   }

   if (node2Ptr == 0) 
   {
      cerr << "NLBeamColumn3d::setDomain: Nd2: ";
      cerr << Nd2 << "does not exist in model\n";
      exit(0);
   }

   // call the DomainComponent class method 
   this->DomainComponent::setDomain(theDomain);
    
   // ensure connected nodes have correct number of dof's
   int dofNode1 = node1Ptr->getNumberDOF();
   int dofNode2 = node2Ptr->getNumberDOF();
   
   if ((dofNode1 !=NND ) || (dofNode2 != NND))
   {
      cerr << "beam2d05::setDomain(): Nd2 or Nd1 incorrect dof ";
      exit(0);
   }
   

   // initialize the transformation
   if (crdTransf->initialize(node1Ptr, node2Ptr))
   {
      cerr << "NLBeamColumn3d::setDomain(): Error initializing coordinate transformation";  
      exit(0);
   }
    
   // get element length
   L = crdTransf->getInitialLength();
   if (L == 0)
   {
      cerr << "NLBeamColumn3d::setDomain(): Zero element length:" << this->getTag();  
      exit(0);
   }
   this->initializeSectionHistoryVariables();
   this->setSectionInterpolation();
}


int
NLBeamColumn3d::commitState()
{
   int err = 0;
   int i = 0;

   do
   {
      vscommit[i] = vs[i];
      err = sections[i++]->commitState();  
   } while (err == 0 && i < nSections);
   
   if (err)
      return err;
   
   // commit the transformation between coord. systems
   if (err = crdTransf->commitState())
      return err;
      
   // commit the element variables state

   distrLoadcommit = prevDistrLoad;
   Uecommit = Uepr;
   kvcommit = kv;
   Secommit = Se;

   return err;
}


int 
NLBeamColumn3d::revertToLastCommit()
{
   int err;
   int i = 0;
      
   do
   {
      vs[i] = vscommit[i];
      err = sections[i]->revertToLastCommit();
      
      Ssr[i] = sections[i]->getStressResultant();
      fs[i]  = sections[i]->getSectionFlexibility();
      
      i++;
   } while (err == 0 && i < nSections);
   
       
   if (err)
      return err;
   
   // revert the transformation to last commit
   if (err = crdTransf->revertToLastCommit())
      return err;
     
   // revert the element state to last commit
   prevDistrLoad = distrLoadcommit;
   Uepr = Uecommit;
   Se   = Secommit;
   kv   = kvcommit;
   
   // compute global resisting forces and tangent
   static Vector currDistrLoad(NL);
   currDistrLoad.Zero();  // SPECIFY LOAD HERE!!!!!!!!! 
   P = crdTransf->getGlobalResistingForce(Se, currDistrLoad);
   K = crdTransf->getGlobalStiffMatrix(kv, Se);
   
   return err;
}


int 
NLBeamColumn3d::revertToStart()
{
   // revert the sections state to start
   int err;
   int i = 0;
     
   do
   {
      fs[i].Zero();
      vs[i].Zero();
      Ssr[i].Zero();
      err = sections[i++]->revertToStart();
 
   } while (err == 0 && i < nSections);

   if (err)
      return err;
   
   // revert the transformation to start
   if (err = crdTransf->revertToStart())
      return err;
  
   // revert the element state to start
   prevDistrLoad.Zero();
   Uepr.Zero();
   Se.Zero();
   kv.Zero();

   P.Zero();
   K.Zero();
   
   return err;
}



const Matrix &
NLBeamColumn3d::getTangentStiff(void)
{
   this->updateElementState();
   return K;
}
    

const Vector &
NLBeamColumn3d::getResistingForce(void)
{
   this->updateElementState();
   return P;
}


void
NLBeamColumn3d::initializeSectionHistoryVariables (void)
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



void
NLBeamColumn3d::setSectionInterpolation (void)
{
    GaussLobattoQuadRule1d01 quadrat(nSections);
    Matrix xi_pt = quadrat.getIntegrPointCoords();   

    for (int i = 0; i < nSections; i++)
    {
	int order = sections[i]->getOrder();
	const ID &code = sections[i]->getType();
	
	b[i] = Matrix(order,NEBD);
	this->getForceInterpolatMatrix(xi_pt(i,0), b[i], code);

	bp[i] = Matrix(order,NL);
	this->getDistrLoadInterpolatMatrix(xi_pt(i,0), bp[i], code);
    }
}



int NLBeamColumn3d::updateElementState(void)
{
  // get element global end displacements
  static Vector Ue(NEGD);
  this->getGlobalDispls(Ue);
 
  // compute global end displacement increments
  static Vector dUe(NEGD);
  // dUe = Ue - Uepr
  dUe = Ue;
  dUe.addVector(1.0, Uepr,-1.0);
  
  if (dUe.Norm() != 0.0  || initialFlag == 0) 
  {
    // compute distributed loads and increments
    static Vector currDistrLoad(NL);
    static Vector distrLoadIncr(NL); 

    currDistrLoad.Zero();  // SPECIFY LOAD HERE!!!!!!!!! 
    distrLoadIncr = currDistrLoad - prevDistrLoad;
    prevDistrLoad = currDistrLoad;

    // update the end displacements
    Uepr = Ue;

    // update the transformation
    crdTransf->update();
       
    // get basic displacements and increments
    static Vector v(NEBD);
    static Vector dv(NEBD);
     
    v = crdTransf->getBasicTrialDisp();    
    dv = crdTransf->getBasicIncrDeltaDisp();    

    // get integration point positions and weights
    int nIntegrPts = nSections;
    // GaussQuadRule1d01 quadrat(nIntegrPts);
    GaussLobattoQuadRule1d01 quadrat(nIntegrPts);
    // const Matrix &xi_pt = quadrat.getIntegrPointCoords();
    const Vector &weight = quadrat.getIntegrPointWeights();
     
    // numerical integration 
    static Vector dSs;       // section internal force increments
    static Vector Ss;        // section "applied" forces (in equilibrium with end forces)
    static Vector dvs;       // section residual deformations
    
    static Vector vr(NEBD);       // element residual displacements
    static Matrix f(NEBD,NEBD);   // element flexibility matrix

    static Matrix I(NEBD,NEBD);   // an identity matrix for matrix inverse
    double dW;                    // section strain energy (work) norm 
    int i;
    
    I.Zero();
    for (i=0; i<NEBD; i++)
      I(i,i) = 1.0;

    // calculate nodal force increments and update nodal forces
    static Vector dSe(NEBD);

    // dSe = kv * dv;
    dSe.addMatrixVector(0.0, kv, dv, 1.0);

    for (int j=0; j < maxIters; j++)
    {
      Se += dSe;
  
      // initialize f and vr for integration
      f.Zero();
      vr.Zero();

      for (i=0; i<nIntegrPts; i++)
      {
	  
          // initialize vectors with correct size  - CHANGE LATER
	  Ss  = Ssr[i];
	  dSs = vs[i];
          dvs = vs[i]; 	  
	  
	  // calculate total section forces
	  // Ss = b*Se + bp*currDistrLoad;
	  Ss.addMatrixVector(0.0, b[i], Se, 1.0);
	  Ss.addMatrixVector(1.0, bp[i], currDistrLoad, 1.0);

	  dSs = Ss;
	  dSs.addVector(1.0, Ssr[i], -1.0);  // dSs = Ss - Ssr[i];
 
	  // compute section deformation increments and update current deformations
	  //       vs += fs * dSs;     
	  
	  dvs.addMatrixVector(0.0, fs[i], dSs, 1.0);
	  
	  vs[i] += dvs;
	  
	  sections[i]->setTrialSectionDeformation(vs[i]);
	  
	  // get section resisting forces
	  Ssr[i] = sections[i]->getStressResultant();
	  
	  // get section flexibility matrix
          fs[i] = sections[i]->getSectionFlexibility();
	  
	  // calculate section residual deformations
	  // dvs = fs * (Ss - Ssr);
	  
	  dSs = Ss;
	  dSs.addVector(1.0, Ssr[i], -1.0);  // dSs = Ss - Ssr[i];
            
	  dvs.addMatrixVector(0.0, fs[i], dSs, 1.0);
	      
          // integrate element flexibility matrix
	  // f = f + (b^ fs * b) * weight(i);
      	  f.addMatrixTripleProduct(1.0, b[i], fs[i], weight(i));

	  // integrate residual deformations
	  // vr += (b^ (vs + dvs)) * weight(i);
	  vr.addMatrixTransposeVector(1.0, b[i], vs[i] + dvs, weight(i));
      }  
      
      f  *= L;
      vr *= L;

      // calculate element stiffness matrix
      // invertMatrix(5, f, kv);

      if (f.Solve(I,kv) < 0)
	 g3ErrorHandler->warning("NLBeamColumn3d::updateElementState() - could not invert flexibility\n");

      dv = v - vr;
      
      // dSe = kv * dv;
      dSe.addMatrixVector(0.0, kv, dv, 1.0);
      
      dW = dv^ dSe;
      if (dW < tol)
        break;
    }     
      
    // determine resisting forces
    Se += dSe;

    // get resisting forces and stiffness matrix in global coordinates
    P = crdTransf->getGlobalResistingForce(Se, currDistrLoad);
    K = crdTransf->getGlobalStiffMatrix(kv, Se);

    initialFlag = 1;

    return 1;
  }
  
  else
    return 0;
}



void NLBeamColumn3d::getGlobalDispls(Vector &dg) const
{
   // determine global displacements
   const Vector &disp1 = node1Ptr->getTrialDisp();
   const Vector &disp2 = node2Ptr->getTrialDisp();

   for (int i = 0; i < NND; i++)
   {
      dg(i)     = disp1(i);
      dg(i+NND) = disp2(i);
   }
}



void NLBeamColumn3d::getGlobalAccels(Vector &ag) const
{
   // determine global displacements
   const Vector &accel1 = node1Ptr->getTrialAccel();
   const Vector &accel2 = node2Ptr->getTrialAccel();

   for (int i = 0; i < NND; i++)
   {
      ag(i)     = accel1(i);
      ag(i+NND) = accel2(i);
   }
}


void NLBeamColumn3d::getForceInterpolatMatrix(double xi, Matrix &b, const ID &code)
{
   b.Zero();
    
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
	    b(i,1) = 1.0/L;
	    b(i,2) = 1.0/L;
	    break;
	 case SECTION_RESPONSE_MY:		// Moment, My, interpolation
	    b(i,3) = xi - 1.0;
	    b(i,4) = xi;
	    break;
	 case SECTION_RESPONSE_VZ:		// Shear, Vz, interpolation
	    b(i,3) = 1.0/L;
	    b(i,4) = 1.0/L;
	    break;
	 case SECTION_RESPONSE_T:		// Torque, T, interpolation
	    b(i,5) = 1.0;
	    break;
	 default:
	    break;
      }	
   }
}


void NLBeamColumn3d::getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code)
{
   bp.Zero();

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
	 case SECTION_RESPONSE_MY:		// Moment, My, interpolation
	    bp(i,2) = xi*(1-xi)*L*L/2;
	    break;
	 case SECTION_RESPONSE_VZ:		// Shear, Vz, interpolation
	    bp(i,2) = (0.5-xi)*L;
	    break;
	 case SECTION_RESPONSE_T:		// Torsion, T, interpolation
	    break;
	 default:
	    break;
      }
   }
}


const Matrix &
NLBeamColumn3d::getSecantStiff(void)
{
    return this->getTangentStiff();
}

    
const Matrix &
NLBeamColumn3d::getDamp(void)
{
    return d;  // zero matrix still 
}


const Matrix &
NLBeamColumn3d::getMass(void)
{ 
    m(0,0) = m(1,1) = m(2,2) =
    m(6,6) = m(7,7) = m(8,8) = rho * L/2;
  
    return m;
}


void 
NLBeamColumn3d::zeroLoad(void)
{
    load.Zero();
}

int
NLBeamColumn3d::addLoad(const Vector &moreLoad)
{
    if (moreLoad.Size() != NEGD) 
    {
	cerr << "NLBeamColumn3d::addLoad: vector not of correct size\n";
	return -1;
    }
    load += moreLoad;
    return 0;
}


const Vector &
NLBeamColumn3d::getResistingForceIncInertia()
{	
    static Vector f(NEGD);
    static Vector ag(NEGD);
    
    f = this->getResistingForce();
    this->getMass();
    this->getGlobalAccels(ag);
    
    // Pinert = f + m * ag;
    Pinert = f;
    Pinert.addMatrixVector(1.0, m, ag, 1.0);
   
    return Pinert;
}




bool
NLBeamColumn3d::isSubdomain(void)
{
    return false;
}

int
NLBeamColumn3d::sendSelf(int commitTag, Channel &theChannel)
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
      if (crdTransfDbTag  != 0) {
	   crdTransf->setDbTag(crdTransfDbTag);
       }
  }
  idData(7) = crdTransfDbTag;
  

  if (theChannel.sendID(dbTag, commitTag, idData) < 0)  
  {
     g3ErrorHandler->warning("NLBeamColumn3d::sendSelf() - %s\n",
	     		     "failed to send ID data");
     return -1;
  }    
  
  if (crdTransf->sendSelf(commitTag, theChannel) < 0)  
  {
     g3ErrorHandler->warning("NLBeamColumn3d::sendSelf() - %s\n",
	     		     "failed to send crdTranf");
     return -1;
  }      

  
  //
  // send an ID for the sections containing each sections dbTag and classTag
  // if section ha no dbTag get one and assign it
  //

  ID idSections(2*nSections);
  loc = 0;
  for (i = 0; i<nSections; i++) 
  {
    int sectClassTag = sections[i]->getClassTag();
    int sectDbTag = sections[i]->getDbTag();
    if (sectDbTag == 0) 
    {
      sectDbTag = theChannel.getDbTag();
      sections[i]->setDbTag(sectDbTag);
    }

    idSections(loc) = sectClassTag;
    idSections(loc+1) = sectDbTag;
    loc += 2;
  }

  if (theChannel.sendID(dbTag, commitTag, idSections) < 0)  {
    g3ErrorHandler->warning("NLBeamColumn3d::sendSelf() - %s\n",
			    "failed to send ID data");
    return -1;
  }    

  //
  // send the sections
  //
  
  for (j = 0; j<nSections; j++) {
    if (sections[j]->sendSelf(commitTag, theChannel) < 0) {
      g3ErrorHandler->warning("NLBeamColumn3d::sendSelf() - section %d %s\n",
			      j,"failed to send itself");
      return -1;
    }
  }
  
  // into a vector place distrLoadCommit, rho, UeCommit, Secommit and kvcommit
  int secDefSize = 0;
  for (i = 0; i < nSections; i++)
  {
     int size = sections[i]->getOrder();
     secDefSize   += size;
  }

  
  
  static Vector dData(1+1+NL+NEGD+NEBD+NEBD*NEBD+secDefSize); 
  loc = 0;

  // place double variables into Vector
  dData(loc++) = rho;
  dData(loc++) = tol;
  
  // put  distrLoadCommit into the Vector
  for (i=0; i<NL; i++) 
  {
    dData(loc) = distrLoadcommit(i);
    loc++;
  }

  // place UeCommit into Vector
  for (i=0; i<NEGD; i++)
    dData(loc++) = Uecommit(i);

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

  
  if (theChannel.sendVector(dbTag, commitTag, dData) < 0)  
  {
     g3ErrorHandler->warning("NLBeamColumn3d::sendSelf() - %s\n",
	 		     "failed to send Vector data");
     return -1;
  }    


  return 0;
}


int
NLBeamColumn3d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
   //
  // into an ID of size 7 place the integer data
  //
  int dbTag = this->getDbTag();
  int i,j,k;
  

  static ID idData(9); // one bigger than needed 

  if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
    g3ErrorHandler->warning("NLBeamColumn3d::recvSelf() - %s\n",
			    "failed to recv ID data");
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
      crdTransf = theBroker.getNewCrdTransf3d(crdTransfClassTag);
      if (crdTransf == 0) {
	  g3ErrorHandler->warning("NLBeamColumn3d::recvSelf() - %s %d\n",
				  "failed to obtain a CrdTrans object with classTag",
				  crdTransfClassTag);
	  return -2;	  
      }
  }
  crdTransf->setDbTag(crdTransfDbTag);
  
  // invoke recvSelf on the crdTransf obkject
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0)  
  {
     g3ErrorHandler->warning("NLBeamColumn3d::sendSelf() - %s\n",
	     		     "failed to recv crdTranf");
     return -3;
  }      

  
  
  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*idData(3));
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    g3ErrorHandler->warning("NLBeamColumn3d::recvSelf() - %s\n",
			    "failed to recv ID data");
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

    // create a new array to hold pointers
    sections = new SectionForceDeformation *[idData(3)];
    if (sections == 0) {
      g3ErrorHandler->fatal("NLBeamColumn3d::recvSelf() - %s %d\n",
			      "out of memory creating sections array of size",idData(3));
      return -1;
    }    

    // create a section and recvSelf on it
    nSections = idData(3);
    loc = 0;
    
    for (i=0; i<nSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;
      sections[i] = theBroker.getNewSection(sectClassTag);
      if (sections[i] == 0) {
	g3ErrorHandler->fatal("NLBeamColumn3d::recvSelf() - %s %d\n",
			      "Broker could not create Section of class type",sectClassTag);
	return -1;
      }
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	g3ErrorHandler->warning("NLBeamColumn3d::recvSelf() - section %d %s\n",
				i,"failed to recv itself");
	return -1;
      }     
    }
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
	  g3ErrorHandler->fatal("NLBeamColumn3d::recvSelf() - %s %d\n",
				"Broker could not create Section of class type",sectClassTag);
	  return -1;
	}
      }

      // recvvSelf on it
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	g3ErrorHandler->warning("NLBeamColumn3d::recvSelf() - section %d %s\n",
				i,"failed to recv itself");
	return -1;
      }     
    }
  }
  
  
  
  // into a vector place distrLoadCommit, rho, UeCommit, Secommit and kvcommit
  int secDefSize = 0;
  for (int ii = 0; ii < nSections; ii++)
  {
     int size = sections[ii]->getOrder();
     secDefSize   += size;
  }
  
  Vector dData(1+1+NL+NEGD+NEBD+NEBD*NEBD+secDefSize);   
  
  if (theChannel.recvVector(dbTag, commitTag, dData) < 0)  {
    g3ErrorHandler->warning("NLBeamColumn3d::sendSelf() - %s\n",
			    "failed to send Vector data");
    return -1;
  }    
  
  loc = 0;
  
  // place double variables into Vector
  rho = dData(loc++);
  tol = dData(loc++);
  
  // put  distrLoadCommit into the Vector
  for (i=0; i<NL; i++) 
  {
    dData(loc) = distrLoadcommit(i) = dData(loc++);
    loc++;
  }

  // place UeCommit into Vector
  for (i=0; i<NEGD; i++)
    Uecommit(i) = dData(loc++);

  // place kvcommit into vector
  for (i=0; i<NEBD; i++) 
    Secommit(i) = dData(loc++);

  // place kvcommit into vector
  for (i=0; i<NEBD; i++) 
     for (j=0; j<NEBD; j++)
        kvcommit(i,j) = dData(loc++);

	prevDistrLoad = distrLoadcommit;
	Uepr = Uecommit;
	kv   = kvcommit;
	Se   = Secommit;

	// Delete the old
	if (vscommit != 0)
		delete [] vscommit;

	// Allocate the right number
	vscommit = new Vector[nSections];
	if (vscommit == 0) {
		g3ErrorHandler->warning("%s -- failed to allocate vscommit array",
			"NLBeamColumn3d::recvSelf");
		return -1;
	}

	for (k = 0; k < nSections; k++) {
		int order = sections[k]->getOrder();

		// place vscommit into vector
		vscommit[k] = Vector(order);
		for (i = 0; i < order; i++)
			(vscommit[k])(i) = dData(loc++);
	}
	
	// Delete the old
	if (fs != 0)
		delete [] fs;

	// Allocate the right number
	fs = new Matrix[nSections];  
	if (fs == 0) {
		g3ErrorHandler->warning("%s -- failed to allocate fs array",
			"NLBeamColumn3d::recvSelf");
		return -1;
	}
   
	// Delete the old
	if (vs != 0)
		delete [] vs;

	// Allocate the right number
	vs = new Vector[nSections];  
	if (vs == 0) {
		g3ErrorHandler->warning("%s -- failed to allocate vs array",
			"NLBeamColumn3d::recvSelf");
		return -1;
	}
	
	// Delete the old
	if (Ssr != 0)
		delete [] Ssr;

	// Allocate the right number
	Ssr = new Vector[nSections];  
	if (Ssr == 0) {
		g3ErrorHandler->warning("%s -- failed to allocate Ssr array",
			"NLBeamColumn3d::recvSelf");
		return -1;
	}

	// Set up section history variables 
	this->initializeSectionHistoryVariables();

	// Delete the old
	if (b != 0)
		delete [] b;

	// Allocate the right number
	b = new Matrix[nSections];  
	if (b == 0) {
		g3ErrorHandler->warning("%s -- failed to allocate b array",
			"NLBeamColumn3d::recvSelf");
		return -1;
	}

	// Delete the old
	if (bp != 0)
		delete [] bp;

	// Allocate the right number
	bp = new Matrix[nSections];  
	if (bp == 0) {
		g3ErrorHandler->warning("%s -- failed to allocate bp array",
			"NLBeamColumn3d::recvSelf");
		return -1;
	}

	// Set up section force interpolation matrices
	this->setSectionInterpolation();

	return 0;
}


void 
NLBeamColumn3d::compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const
{
   // update the transformation
   crdTransf->update();
       
   // get basic displacements and increments
   static Vector ub(NEBD);
   ub = crdTransf->getBasicTrialDisp();    
  
   // get integration point positions and weights
   int nIntegrPts = nSections;
   GaussLobattoQuadRule1d01 quadrat(nIntegrPts);
   const Matrix &xi_pt  = quadrat.getIntegrPointCoords();

   // setup Vandermode and CBDI influence matrices
   int i;
   double xi;
 
   // get CBDI influence matrix
   Matrix ls(nIntegrPts, nIntegrPts);
   getCBDIinfluenceMatrix(nIntegrPts, xi_pt, L, ls);
     
   // get section curvatures
   Vector kappa_y(nIntegrPts);  // curvature
   Vector kappa_z(nIntegrPts);  // curvature
   static Vector vs;                // section deformations 
        
   for (i=0; i<nSections; i++)
   {
       // THIS IS VERY INEFFICIENT ... CAN CHANGE IF RUNS TOO SLOW
       int sectionKey1 = 0;
       int sectionKey2 = 0;
       const ID &code = sections[i]->getType();
       int j;
       for (j = 0; j < code.Size(); j++)
       {
	   if (code(j) == SECTION_RESPONSE_MZ)
	       sectionKey1 = j;
	   if (code(j) == SECTION_RESPONSE_MY)
	       sectionKey2 = j;
       }
       if (sectionKey1 == 0)
	   g3ErrorHandler->fatal("FATAL NLBeamColumn3d::compSectionResponse - section does not provide Mz response\n");
       if (sectionKey2 == 0)
	   g3ErrorHandler->fatal("FATAL NLBeamColumn3d::compSectionResponse - section does not provide My response\n");

       // get section deformations
       vs = sections[i]->getSectionDeformation();
       
       kappa_z(i) = vs(sectionKey1);
       kappa_y(i) = vs(sectionKey2); 
   }

   //cout << "kappa_y: " << kappa_y;   
   //cout << "kappa_z: " << kappa_z;   
   
   Vector v(nIntegrPts), w(nIntegrPts);
   static Vector xl(NDM), uxb(NDM);
   static Vector xg(NDM), uxg(NDM); 
   // double theta;                             // angle of twist of the sections

   // v = ls * kappa_z;  
   v.addMatrixVector (0.0, ls, kappa_z, 1.0);  
   // w = ls * kappa_y *  (-1);  
   w.addMatrixVector (0.0, ls, kappa_y, -1.0);
   
   for (i=0; i<nSections; i++)
   {
      xi = xi_pt(i,0);

      xl(0) = xi * L;
      xl(1) = 0;
      xl(2) = 0;

      // get section global coordinates
      sectionCoords[i] = crdTransf->getPointGlobalCoordFromLocal(xl);

      // compute section displacements
      //theta  = xi * ub(5); // consider linear variation for angle of twist. CHANGE LATER!!!!!!!!!!
      uxb(0) = xi * ub(0); // consider linear variation for axial displacement. CHANGE LATER!!!!!!!!!!
      uxb(1) = v(i);
      uxb(2) = w(i);
          
      // get section displacements in global system 
      sectionDispls[i] = crdTransf->getPointGlobalDisplFromBasic(xi, uxb);
   }	       
}



void
NLBeamColumn3d::Print(ostream &s, int flag)
{
   if (flag == 1)
   {
/*    static Vector yAxis(3), zAxis(3);
     
      this->getLocalAxes(yAxis, zAxis);
                        
      s << "\n#ELEMENT " << xAxis(0) << " " << xAxis(1) << " " << xAxis(2)
                  << " " << zAxis(0) << " " << zAxis(1) << " " << zAxis(2); 
   
      // print node coordinates
      const Vector &node1Crd = node1Ptr->getCrds();
      const Vector &node2Crd = node2Ptr->getCrds();

      // determine global displacements
      const Vector &node1Disp = node1Ptr->getTrialDisp();
      const Vector &node2Disp = node2Ptr->getTrialDisp();

      s << "\n#NODE " << node1Crd(0) << " " << node1Crd(1) << " " << node1Crd(2)
               << " " << node1Disp(0) << " " << node1Disp(1) << " " << node1Disp(2)
               << " " << node1Disp(3) << " " << node1Disp(4) << " " << node1Disp(5);

      s << "\n#NODE " << node2Crd(0) << " " << node2Crd(1) << " " << node2Crd(2)
               << " " << node2Disp(0) << " " << node2Disp(1) << " " << node2Disp(2)
               << " " << node2Disp(3) << " " << node2Disp(4) << " " << node2Disp(5);
      
      this->compSectionResponse();

      static Vector S(2);  // shear forces and torsion;
      static Vector currDistrLoad(NL);
  
      currDistrLoad.Zero();  // SPECIFY LOAD HERE!!!!!!!!! 

      for (i=0; i<nSections; i++)
      {
         xi = xi_pt(i,0);
 
         S(0) =  (Se(2)+Se(1))/L + currDistrLoad(1)*(xi - 0.5)* L;
         S(1) = -(Se(3)+Se(4))/L + currDistrLoad(2)*(xi - 0.5)* L;
 
         sections[i]->setShearForces(S);
      } 

      
      for (int i = 0; i < nSections; i++)
         sections[i]->Print(s, flag); */
   }
   else
   {
      s << "\nElement: " << this->getTag() << " Type: NLBeamColumn3d ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << nSections;

//    for (int i = 0; i < nSections; i++)
//       s << "\nSection "<<i<<" :" << *sections[i];

//       s << "\nSection "<<0<<" :" << *sections[0];

       s << "\tStiffness Matrix:\n" << kv;
         s << "\tResisting Force: " << Se;

   }
}


ostream &operator<<(ostream &s, NLBeamColumn3d &E)
{
    E.Print(s);
    return s;
}

int
NLBeamColumn3d::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
   
   if (displayMode == 1) 
   {
       // first determine the two end points of the element based on
       //  the display factor (a measure of the distorted image)
    
       static Vector v1(NDM), v2(NDM);

       const Vector &node1Crd = node1Ptr->getCrds();
       const Vector &node2Crd = node2Ptr->getCrds();	
       const Vector &node1Disp = node1Ptr->getDisp();
       const Vector &node2Disp = node2Ptr->getDisp();    

       int i;
       
       // allocate array of vectors to store section coordinates and displacements
       
       Vector *coords = new Vector [nSections];
     
       if (!coords)
       {
	   cerr << "NLBeamColumn3d::displaySelf() -- failed to allocate coords array";   
	   exit(-1);
       }
       
       for (i = 0; i < nSections; i++)
	  coords[i] = Vector(NDM);
       
       Vector *displs = new Vector [nSections];
     
       if (!displs)
       {
	   cerr << "NLBeamColumn3d::displaySelf() -- failed to allocate coords array";   
	   exit(-1);
       }

       for (i = 0; i < nSections; i++)
	  displs[i] = Vector(NDM);
       
       int error;
       
       this->compSectionDisplacements(coords, displs);

       // get global displacements and coordinates of each section          

       v1 = node1Crd + node1Disp*fact;
 
       // get global displacements and coordinates of each section          

       for (i = 0; i<nSections; i++) 
       {
          v2 = coords[i] + displs[i]*fact;
       
          error = theViewer.drawLine(v1, v2, 1.0, 1.0);

          if (error)
            return error;
          v1 = v2;

       }  
       
       v2 = node2Crd + node2Disp*fact;
       
       error = theViewer.drawLine(v1, v2, 1.0, 1.0);

       if (error)
	  return error;
   }
   return 0;
}




int 
NLBeamColumn3d::setResponse(char **argv, int argc, Information &eleInformation)
{
    //
    // we compare argv[0] for known response types 
    //

    // force - 
    if ((strcmp(argv[0],"forces") == 0) || (strcmp(argv[0],"force") == 0)) {
	Vector *newVector = new Vector(NEBD);
	if (newVector == 0) {
	    cerr << "NLBeamColumn3d::setResponse() - out of memory creating vector\n";
	    return -1;
	}	
	eleInformation.theVector = newVector;	
	eleInformation.theType = VectorType;
	return 1;
    } 

    // section response -
    else if (strcmp(argv[0],"section") ==0) {
	if (argc <= 2)
	    return -1;
	
	int sectionNum = atoi(argv[1]);
	if ((sectionNum > 0) && (sectionNum <= nSections)) {
	    int ok = sections[sectionNum-1]->setResponse(&argv[2], argc-2, 
							 eleInformation);
	    if (ok < 0)
		return -1;
	    else if (ok >= 0 && ok < MAX_SECTION_RESPONSE_ID)
		return sectionNum*MAX_SECTION_RESPONSE_ID + ok;  
	    else 
		return -1;
	}
    }
    
    return -1;
}

int 
NLBeamColumn3d::getResponse(int responseID, Information &eleInformation)
{
  switch (responseID) {

    case -1: // unknown 
      return -1;
      
    case 1:  // forces
      if (eleInformation.theVector != 0)
	  *(eleInformation.theVector) = Se;
      return 0;            

    default: 
      if (responseID >= MAX_SECTION_RESPONSE_ID) { // section quantity
	  int sectionNum = responseID/MAX_SECTION_RESPONSE_ID; 
	  if ((sectionNum > 0) && (sectionNum <= nSections)) {
	      return sections[sectionNum-1]->
		  getResponse(responseID-MAX_SECTION_RESPONSE_ID*sectionNum, 
			      eleInformation);
	  } else
	      return -1;
      } else // unknown
	  return -1;
  }
}


