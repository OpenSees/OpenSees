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

// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLengthInterface2D.cpp,v $
// $Revision: 1.1 $

// Written: Roozbeh Geraili Mikola  (roozbehg@berkeley.edu)
//          Prof. Nicholas Sitar (nsitar@ce.berkeley.edu)
//
// Date: July 02 2010

/*
 - element zeroLengthInterface2D eleTag? -sNdNum sNdNum? -pNdNum pNdNum? –dof sdof? mdof? -Nodes Nodes? Kn? Kt? phi?
 Description: This file contains the implementation for the ZeroLengthInterface2D class.
*/

#include "ZeroLengthInterface2D.h"
#include <Information.h>
#define PI 3.141592653589793238462643383279502884197169399
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>


//*********************************************************************
//  Full Constructor:

ZeroLengthInterface2D::ZeroLengthInterface2D(int tag, int sNdNum, int pNdNum, int sDof, int pDof, const ID& Nodes,
					     double Knormal, double Ktangent, double frictionAngle)
  :Element(tag,ELE_TAG_ZeroLengthInterface2D),
   connectedExternalNodes(sNdNum + pNdNum),
   N(6), T(6), ContactNormal(2), Ki(0), load(0)
{
  //static data
  SecondaryNodeNum = sNdNum;
  PrimaryNodeNum = pNdNum;
  numberNodes = SecondaryNodeNum + PrimaryNodeNum;
  SecondaryDof = sDof;
  PrimaryDof = pDof;
  
  // allocate contact force vector
  pressure.resize(numberNodes);
  // allocate normal gap vector
  normal_gap.resize(numberNodes);
  // allocate shear gap vector
  shear_gap.resize(numberNodes);
  stored_shear_gap.resize(numberNodes);
  
  // set the vectors to zero
  for(int i = 0; i < numberNodes; i++)
    {
      stored_shear_gap(i) = 0;
      shear_gap(i) = 0;
      pressure(i) = 0;
      normal_gap(i) = 0;
    }
  
  int alloc = SecondaryDof * SecondaryNodeNum + PrimaryDof * PrimaryNodeNum;
  
  // static data for 2D
  stiff.resize(alloc, alloc);
  resid.resize(alloc);
  zeroMatrix.resize(alloc, alloc);
  
  /*
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 3)
  opserr << "FATAL ZeroLength::setUp - failed to create an ID of correct size\n";
  */
  
  // create an ID of correct size
  nodePointers = new Node* [numberNodes];
  
  // restore node number
  for (int i = 0 ; i < numberNodes; i++) connectedExternalNodes(i) = (int) Nodes(i);
  
  // assign Kn, Kt, fc
  Kn = Knormal;
  Kt = Ktangent;
  // friction ration fc = tan(phi0
  fc = tan (frictionAngle * PI / 180);
  
  // initialized contact flag be zero
  ContactFlag=0;
  
}

//null constructor

ZeroLengthInterface2D::ZeroLengthInterface2D(void)                               //fixme numberNodes?
  :Element(0,ELE_TAG_ZeroLengthInterface2D),
  connectedExternalNodes(numberNodes),
  N(6), T(6), Ki(0), load(0)
{
    // ensure the connectedExternalNode ID is of correct size
    if (connectedExternalNodes.Size() != numberNodes)
      opserr << "FATAL ZeroLengthInterface2D::ZeroLengthInterface2D - failed to create an ID of correct size\n";
    for (int j = 0; j < numberNodes; j++ ) nodePointers[j] = 0;
}


//  Destructor:
//  delete must be invoked on any objects created by the object
//  and on the matertial object.
ZeroLengthInterface2D::~ZeroLengthInterface2D()
{
  if (load != 0) delete load;
  if (Ki != 0) delete Ki;
}

int
ZeroLengthInterface2D::getNumExternalNodes(void) const
{
    return numberNodes;
}

const ID &
ZeroLengthInterface2D::getExternalNodes(void)
{
  return connectedExternalNodes;
}

Node **
ZeroLengthInterface2D::getNodePtrs(void)
{
    return nodePointers;
}

int
ZeroLengthInterface2D::getNumDOF(void)
{
     return numDOF;
}

// method: setDomain()
// to set a link to the enclosing Domain and to set the node pointers.
// also determines the number of dof associated
// with the ZeroLengthInterface2D element
void
ZeroLengthInterface2D::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    for (int j = 0; j < numberNodes; j++ )  nodePointers[j] = 0;
    return;
  }
  
  // call the base class method
  this->DomainComponent::setDomain(theDomain);
  
  numDOF = 0;
  // set default values for error conditions
  // first set the node pointers
  for (int i = 0; i < numberNodes; i++ ) {
    int Nd = connectedExternalNodes(i);
    nodePointers[i] = theDomain->getNode(Nd);
    if ( nodePointers[i] == 0 ) {
      opserr << "WARNING ZeroLengthInterface2D::setDomain() - Nd: " << Nd << " does not exist in ";
      return; 
    }
    // now determine the number of dof and the dimension
    int dofNd = nodePointers[i]->getNumberDOF();
    // if differing dof at the ends - print a warning message
    //if ( dofNd != 2 || dofNd != 3) {
    //	opserr << "WARNING ZeroLengthInterface2D::setDomain cannot handle " << dofNd << "dofs\n";
    //    return;
    //}
    numDOF += dofNd;
  }
  
  /*
  // Check that length is zero within tolerance
  const Vector &end1Crd = nodePointers[0]->getCrds();
  const Vector &end2Crd = nodePointers[1]->getCrds();
  const Vector &end2Crd = nodePointers[2]->getCrds();
  Vector diff = end1Crd - end2Crd;
  double L  = diff.Norm();
  double v1 = end1Crd.Norm();
  double v2 = end2Crd.Norm();
  double vm;
  
  vm = (v1<v2) ? v2 : v1;
  
  if (L > LENTOL*vm)
  opserr << "WARNING ZeroLengthInterface2D::setDomain(): Element " << this->getTag() << " has L= " << L <<
  ", which is greater than the tolerance\n";
  */
}

int
ZeroLengthInterface2D::commitState()
{
  // need to update stick point here
  if (ContactFlag == 2 )   // slide case, update stick point
    for(int i = 0; i < numberNodes; i++) 
      stored_shear_gap(i) = shear_gap(i);
  
  return 0;
}

int
ZeroLengthInterface2D::revertToLastCommit()
{
  ///////////////////////////////////////////
  // need to revert the stickPoint??
  //int code=0;
  // revert state for 1d materials
  //for (int i=0; i<numMaterials1d; i++)
  //code += theMaterial1d[i]->revertToLastCommit();
  //return code;
  
  //xi=stickPt;
  return 0;
  
}

int
ZeroLengthInterface2D::revertToStart()
{
  // need to rezero stickPoint??
  //int code=0;
  // revert to start for 1d materials
  //for (int i=0; i<numMaterials1d; i++)
  //code += theMaterial1d[i]->revertToStart();
  //return code;
  // zero stickPt
  
  for(int i = 0; i < numberNodes; i++) 
    stored_shear_gap(i) = 0;
  
  return 0;
}

// update
// calculate stress-strain relation -- M. Frank
int
ZeroLengthInterface2D::update(void)
{
  return 0;
}

const Matrix &
ZeroLengthInterface2D::getTangentStiff(void)
{
  int tang_flag = 1 ; //get the tangent
  //zero stiffness and residual
  stiff.Zero();
  //do tangent and residual here
  formGlobalResidAndTangent( tang_flag ) ;
  //opserr<< stiff ;
  return stiff ;
}

const Matrix &
ZeroLengthInterface2D::getInitialStiff(void)
{
  int tang_flag = 1 ; //get the tangent
  stiff.Zero();
  //do tangent and residual here
  formGlobalResidAndTangent( tang_flag ) ;
  return stiff ;
}

const Matrix &
ZeroLengthInterface2D::getDamp(void)
{
  // no damp
  zeroMatrix.Zero();
  return zeroMatrix;
}

const Matrix &
ZeroLengthInterface2D::getMass(void)
{
  // no mass
  zeroMatrix.Zero();
  return zeroMatrix;
}

void
ZeroLengthInterface2D::zeroLoad(void)
{
  // does nothing now
}

int
ZeroLengthInterface2D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // meaningless to addLoad to a contact !
  return 0;
}

int
ZeroLengthInterface2D::addInertiaLoadToUnbalance(const Vector &accel)
{
  // does nothing as element has no mass yet!
  return 0;
}

const Vector &
ZeroLengthInterface2D::getResistingForce()
{
  int tang_flag = 0 ; //don't get the tangent
  resid.Zero();
  formGlobalResidAndTangent( tang_flag ) ;
  //opserr<< "resid="<<resid;
  return resid ;
}

const Vector &
ZeroLengthInterface2D::getResistingForceIncInertia()
{
  // there is no Inertia
  int tang_flag = 0 ; //don't get the tangent
  resid.Zero();
  formGlobalResidAndTangent( tang_flag ) ;
  return  resid ;
}

int
ZeroLengthInterface2D::sendSelf(int commitTag, Channel &theChannel)
{
  // doing nothing here
  return 0;
}

int
ZeroLengthInterface2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  // doing nothing here
  return 0;
}

int
ZeroLengthInterface2D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{ 
  // nothing to display
  return 0;
}

void
ZeroLengthInterface2D::Print(OPS_Stream &s, int flag)
{
  if (flag == 0) { // print everything
    s << "Element: " << this->getTag();
    s << " type: ZeroLengthInterface2D  Nodes: " << connectedExternalNodes << endln;
  } else if (flag == 1) {
    s << this->getTag() << "  ";
  }
}

Response*
ZeroLengthInterface2D::setResponse(const char **argv, int argc, Information &eleInformation)
{
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
    return new ElementResponse(this, 1, resid);
  // tangent stiffness matrix
  else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0)
    return new ElementResponse(this, 2, stiff);
  // contact pressure
  else if (strcmp(argv[0],"pressure")== 0)
    {
      return new ElementResponse(this, 3, pressure);
    } else if (strcmp(argv[0],"gap")== 0)
    return new ElementResponse(this, 4, normal_gap);
  else
    return 0;
}

int
ZeroLengthInterface2D::getResponse(int responseID, Information &eleInfo)
{
 if (responseID == 1)
	 return eleInfo.setVector(this->getResistingForce());
 else if (responseID == 2)
	 return eleInfo.setMatrix(this->getTangentStiff());
 else if (responseID == 3)
	 return eleInfo.setVector(this->pressure);
 else if (responseID == 4)
	 return eleInfo.setVector(this->normal_gap);
 else
	 return -1;
}

// Private methods
// determine the secondary/primary pair in contact, and setup Vectors (N,T1,T2)
int ZeroLengthInterface2D::contactDetect(int s, int p1, int p2, int stage)
{
  //+--------------+-----------------+----------------+----------------+---------------+
  // NOTES: some methods to get displacements from nodes
  //+--------------+-----------------+----------------+----------------+---------------+
  // getDisp() :         get commit(k-1) disp, will be commit(k) after commit
  // getTrialDisp():     get Trial(k) disp
  // getIncrDisp():      get Trial(k)-Commit(k-1), will be 0 after commit
  // getIncrDeltaDisp(): get Trial(k)-Trial(k-1),  will be 0 after commit
  //+--------------+-----------------+----------------+----------------+--------------
  ////////////////////////////// for transient gap ///////////////////////////////////
  // DEFINE:
  // gap = (U_primary-U_secondary) / dot(ContactNormal),
  // defines overlapped normal distance, always keep positive (+) when contacted
  ///*
  // get current position and after trial displacement for (secondary, primary1, primary2) nodes
  int i;
  const Vector &xs = nodePointers[s]->getCrds();
  const Vector &uxs = nodePointers[s]->getTrialDisp();
  const Vector &x1 = nodePointers[p1]->getCrds();
  const Vector &ux1= nodePointers[p1]->getTrialDisp();
  const Vector &x2 = nodePointers[p2]->getCrds();
  const Vector &ux2= nodePointers[p2]->getTrialDisp();
  
  Vector trial_secondary(2), trial_primary1(2), trial_primary2(2);
  for (i = 0; i < 2; i++) {
    trial_secondary(i) = xs(i) + uxs(i);
    trial_primary1(i) = x1(i) + ux1(i);
    trial_primary2(i) = x2(i) + ux2(i);
    //opserr << "trial_secondary: " << trial_secondary(i) << "\n";
    //opserr << "trial_primary1: " << trial_primary1(i) << "\n";
    //opserr << "trial_primary2: " << trial_primary2(i) << "\n";
  }
  
  // calculate normal gap for contact
  Vector diff(2);
  Vector ContactTangent(2);
  for (i = 0; i < 2; i++) {
    diff(i) = trial_primary2(i) - trial_primary1(i);
    //opserr << "diff: " << diff(i) << "\n";
  }
  double L  = diff.Norm();
  // tangent vector
  for (i = 0; i < 2; i++) ContactTangent(i) = (1/L) * (trial_primary2(i) - trial_primary1(i));
  // normal vector
  ContactNormal(0) = - ContactTangent(1);
  ContactNormal(1) = ContactTangent(0);
  
  normal_gap(s) = 0;
  double alpha = 0;
  double alpha_bar = 0;
  for (i = 0; i < 2; i++) {
    alpha += (1/L) * (trial_secondary(i) - trial_primary1(i)) * ContactTangent(i);
    normal_gap(s) += (trial_secondary(i) - trial_primary1(i)) * ContactNormal(i);
    diff(i) = x2(i) - x1(i);
  }
  
  double gapgap = normal_gap(s);
  
  double L_bar = diff.Norm();
  for (i = 0; i < 2; i++) alpha_bar += (1/L_bar) * (xs(i) - x1(i)) * ContactTangent(i);
  shear_gap(s) = (alpha - alpha_bar) * L_bar;
  /*
    /////////////////////////////// for transient gap ///////////////////////////////
    // we have another way to define the gap, can replace previous code block if want
    ////////////////////////////// for dynamic gap //////////////////////////////////
    const Vector   // get current trial incremental position
    &U_secondary = nodePointers[0]->getCrds() + nodePointers[0]->getIncrDisp();
    const Vector
    &U_primary= nodePointers[1]->getCrds() + nodePointers[1]->getIncrDisp();
    gap=0;
    int i;
    for (i=0; i<2; i++){
    gap += (U_primary(i)-U_secondary(i))* ContactNormal(i);
    }
    gap+=gap_n;
    ///////////////// for dynamic gap //////////////////////
    */
  // stage = 0 means searching secondary nodes against primary segments
  // stage = 1 means searching primary nodes against secondary segments
  if ((stage == 0  && normal_gap(s) >= 0 && alpha > 0 && alpha < 1) ||
      (stage == 1  && normal_gap(s) >= 0 && alpha >= 0 && alpha <= 1)) { // in contact
    N(0) = ContactNormal(0);
    N(1) = ContactNormal(1);
    N(2) = -(1 - alpha) * N(0);
    N(3) = -(1 - alpha) * N(1);
    N(4) = -(alpha) * N(0);
    N(5) = -(alpha) * N(1);
    
    T(0) = ContactTangent(0);
    T(1) = ContactTangent(1);
    T(2) = -(1-alpha) * T(0);
    T(3) = -(1-alpha) * T(1);
    T(4) = -(alpha) * T(0);
    T(5) = -(alpha) * T(1);
    
    return 1;
  } else {
    return 0; // Not in contact
  }
}

/*
void  ZeroLengthInterface2D::GlobalResidAndTangentOrder2(int secondary, int primary1, int primary2)
{
    // create a vector for converting local matrix to global
	int sdofNd = nodePointers[secondary]->getNumberDOF();
	int p1dofNd = nodePointers[primary1]->getNumberDOF();
	int p2dofNd = nodePointers[primary2]->getNumberDOF();

	if (sdofNd == 2 && p1dofNd == 3 && p2dofNd == 3) {
        loctoglob[0] = (2 * secondary);
		loctoglob[1] = (2 * secondary) + 1;
	    loctoglob[2] = (2 * SecondaryNodeNum) + (3 * (primary1 - SecondaryNodeNum)); 
		loctoglob[3] = (2 * SecondaryNodeNum) + (3 * (primary1 - SecondaryNodeNum)) + 1; 
	    loctoglob[4] = (2 * SecondaryNodeNum) + (3 * (primary2 - SecondaryNodeNum)); 
		loctoglob[5] = (2 * SecondaryNodeNum) + (3 * (primary2 - SecondaryNodeNum)) + 1;
	} else if (sdofNd == 3 && p1dofNd == 2 && p2dofNd == 2) {
        loctoglob[0] = (2 * SecondaryNodeNum) + (3 * (secondary - SecondaryNodeNum));
		loctoglob[1] = (2 * SecondaryNodeNum) + (3 * (secondary - SecondaryNodeNum)) + 1;
	    loctoglob[2] = (2 * primary1); 
		loctoglob[3] = (2 * primary1) + 1; 
	    loctoglob[4] = (2 * primary2); 
		loctoglob[5] = (2 * primary2) + 1;
	} else if (sdofNd == 2 && p1dofNd == 2 && p2dofNd == 2) {
		loctoglob[0] = (2 * secondary); 
		loctoglob[1] = (2 * secondary) + 1; 
	    loctoglob[2] = (2 * primary1); 
		loctoglob[3] = (2 * primary1) + 1; 
	    loctoglob[4] = (2 * primary2); 
		loctoglob[5] = (2 * primary2) + 1; 
	} else if (sdofNd == 3 && p1dofNd == 3 && p2dofNd == 3) {
	    loctoglob[0] = (3 * secondary); 
		loctoglob[1] = (3 * secondary) + 1; 
	    loctoglob[2] = (3 * primary1); 
		loctoglob[3] = (3 * primary1) + 1; 
	    loctoglob[4] = (3 * primary2); 
		loctoglob[5] = (3 * primary2) + 1; 
	}

}
*/
void  ZeroLengthInterface2D::GlobalResidAndTangentOrder(int secondary, int primary1, int primary2)
{
    // create a vector for converting local matrix to global
	int sdofNd = nodePointers[secondary]->getNumberDOF();
	int p1dofNd = nodePointers[primary1]->getNumberDOF();
	int p2dofNd = nodePointers[primary2]->getNumberDOF();
	int nd[3] = { secondary, primary1, primary2 };
	int dof[3] = { sdofNd, p1dofNd, p2dofNd };
	
	for (int i = 0 ; i < 3 ; i++)
	{
		if (dof[i] == SecondaryDof)
		{
			loctoglob[2 * i]     = (dof[i] * nd[i]);
			loctoglob[2 * i + 1] = (dof[i] * nd[i] + 1);
		} else if (dof[i] == PrimaryDof)
		{
			int add = SecondaryDof * SecondaryNodeNum;
			int pos = nd[i] - SecondaryNodeNum;
			loctoglob[2 * i]     = (dof[i] * pos) + add;
			loctoglob[2 * i + 1] = (dof[i] * pos + 1) + add;
		} else
		{
			// error message dof is no equal with either of secondary and primary dofs
		}
	}
}

void  ZeroLengthInterface2D::formLocalResidAndTangent( int tang_flag , int secondary, int primary1, int primary2, int stage)
{
	// trial frictional force vectors (in local coordinate)
    double t_trial;
    double TtrNorm;
    // Coulomb friction law surface
	double Phi;
	int i, j;

	// set the first value to zero
	pressure(secondary) = 0;
    t_trial=0;

	// int IsContact;
	// detect contact and set flag
    ContactFlag = contactDetect(secondary,primary1,primary2, stage);

	if (ContactFlag == 1) // contacted
	{
	    // create a vector for converting local matrix to global
        GlobalResidAndTangentOrder(secondary, primary1, primary2);

		// contact presure;
	    pressure(secondary) = Kn * normal_gap(secondary);  // pressure is positive if in contact

	    double ng = normal_gap(secondary);

		t_trial = Kt * (shear_gap(secondary) - stored_shear_gap(secondary));  // trial shear force

        // Coulomb friction law, trial state
		//TtrNorm=t_trial.Norm();
		TtrNorm = sqrt(t_trial * t_trial);
		Phi = TtrNorm - fc * pressure(secondary);

		if (Phi <= 0 ) { // stick case
 			if ( tang_flag == 1 ) {
		    // stiff
				for (i = 0; i < 6; i++) {
					for (j = 0; j < 6; j++) {
                        stiff(loctoglob[i],loctoglob[j]) +=   Kn * (N(i) * N(j)) + Kt * (T(i) * T(j));	//2D
					}
				}
			} //endif tang_flag
			// force
			for (i = 0; i < 6; i++) resid(loctoglob[i]) += pressure(secondary) * N(i) + t_trial * T(i);    //2D
		} // end if stick
		else {           // slide case, non-symmetric stiff
            ContactFlag=2;  // set the contactFlag for sliding
            if ( tang_flag == 1 ) {
				// stiff
				for (i = 0; i < 6; i++) {
					for (j = 0; j < 6; j++) {
					stiff(loctoglob[i],loctoglob[j]) += Kn * (N(i) * N(j)) - fc * Kn * (t_trial / TtrNorm) * T(i) * N(j); //2D
					} //endfor i
				} //endfor j
   			    // force
			} // endif tang_flag
            double shear = fc * pressure(secondary) * (t_trial/TtrNorm);
			for (i = 0; i < 6; i++) resid(loctoglob[i]) += (pressure(secondary) * N(i)) + (shear * T(i)) ;      //2D
		} //endif slide
	}  // endif ContactFlag==1
}

void  ZeroLengthInterface2D::formGlobalResidAndTangent( int tang_flag )
{
	// in the first loop the node to node contact will not be considered 
	// in order to prevent node-node contact duplication 
	// but on contrary in the second loop the node to node contact 
	// will be considered and this can be controlled by "stage = 0 or 1"

	// loop over secondary nodes and find the nodes 
    // which are in contact with primary's segments
	for (int i = 0 ; i < SecondaryNodeNum; i++) {
		for (int j = SecondaryNodeNum ; j < SecondaryNodeNum + PrimaryNodeNum - 1; j++) {
			formLocalResidAndTangent( tang_flag, i, j, j+1 , 0);  // stage = 0 //
        }
    }

    // loop over primary nodes and find the nodes 
    // which are in contact with secondary's segments
	for (int i = SecondaryNodeNum ; i < SecondaryNodeNum + PrimaryNodeNum; i++) {
		for (int j = 0 ; j < SecondaryNodeNum - 1; j++) {
            formLocalResidAndTangent( tang_flag, i, j, j+1 , 1);  // stage = 1 //
        }
    }
}


// element zeroLengthInterface2D $eleTag -sNdNum $sNdNum -pNdNum $pNdNum -dof $sdof $mdof -Nodes $Nodes $Kn $Kt $phi
//

static int numZeroLengthInterface2D = 0;
#include <elementAPI.h>

void * 
OPS_ZeroLengthInterface2D(void) {

  if (numZeroLengthInterface2D == 0) {
    numZeroLengthInterface2D++;
    opserr << "ZeroLengthContactNTS2d - Written by Roozbeh G. Mikola and N.Sitar, UC Berkeley\n";
  }

  Element *theEle = 0;
  int numData = 0;  

  // get the ele tag
  int eleTag, sNdNum, pNdNum, sDOF, mDOF;
  numData = 1;

  if (OPS_GetInt(&numData, &eleTag) != 0) {
    opserr << "ZeroLengthInterface2D::WARNING invalid eleTag \n";
    return 0;
  }

  const char *nextString = OPS_GetString();

  if (strcmp(nextString,"-sNdNum") != 0) {
    opserr << "ZeroLengthInterface2D:: expecting -sNdNum \n";
    return 0;
  }

  // get the number of secondary nodes
  numData = 1;
  if (OPS_GetInt(&numData, &sNdNum) != 0) {
    opserr << "ZeroLengthInterface2D::WARNING invalied sNdNum \n";
    return 0;
  }

  numData = 10;
  nextString = OPS_GetString();

  if (strcmp(nextString,"-mNdNum") != 0 && strcmp(nextString,"-pNdNum") != 0) {
    opserr << "ZeroLengthInterface2D:: expecting -pNdNum\n";
    return 0;
  }
  
  numData = 1;
  if (OPS_GetInt(&numData, &pNdNum) != 0) {
    opserr << "ZeroLengthInterface2D::WARNING invalied pNdNum \n";
    return 0;
  }

  numData = 10;
  nextString = OPS_GetString();

  if (strcmp(nextString,"-dof") != 0) {
    opserr << "ZeroLengthInterface2D:: expecting -sdof in "<<
      "element zeroLengthInterface2D eleTag? -sNdNum sNdNum? -pNdNum pNdNum? -dof sdof? mdof? -Nodes Nodes? Kn? Kt? phi? \n" ;
    return 0;
  }

  numData = 1;
  if (OPS_GetInt(&numData, &sDOF) != 0) {
    opserr << "ZeroLengthInterface2D::WARNING invalied sDOF\n";
    return 0;
  }

  numData = 1;
  if (OPS_GetInt(&numData, &mDOF) != 0) {
    opserr << "ZeroLengthInterface2D::WARNING invalied mDOF\n";
    return 0;
  }

  // a quick check on number of args
  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 3 + sNdNum + pNdNum) {
    opserr << "ZeroLengthInterface2D::WARNING too few arguments " <<
      "element zeroLengthInterface2D eleTag? -sNdNum sNdNum? -pNdNum pNdNum? -dof sdof? mdof? -Nodes Nodes? Kn? Kt? phi? \n" ;
    return 0;
  }

  numData = 10;
  nextString = OPS_GetString();

  if (strcmp(nextString,"-Nodes") != 0) {
    opserr << "ZeroLengthInterface2D:: expecting -Nodes\n";
    return 0;
  }

  // read the Nodes values
  numData = sNdNum+pNdNum;
  int *theNodeData = new int[numData];
  ID  Nodes(theNodeData, numData);

  if (OPS_GetInt(&numData, theNodeData) != 0) {
    opserr << "ZeroLengthInterface2D:: not enough node tags provided for ele: ";
    opserr << eleTag << "\n";
    return 0;
  }

  // read the material properties
  numData = 3;
  double dData[3];
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "ZeroLengthInterface2D::WARNING invalid Kn,Kt or phi\n" ;
    return 0;
  }

  //
  // now we create the element and add it to the domain
  //
  
  theEle = new ZeroLengthInterface2D(eleTag, sNdNum, pNdNum, sDOF, mDOF, Nodes, dData[0], dData[1], dData[2]);
  return theEle;
}


