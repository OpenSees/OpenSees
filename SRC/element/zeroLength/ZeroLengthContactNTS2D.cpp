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

// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLengthContactNTS2D.cpp,v $
// $Revision: 1.1 $

// Written: Roozbeh Geraili Mikola  (roozbehg@berkeley.edu)
//          Nicholas Sitar          (sitar@ce.berkeley.edu)
//
// Created: 12/28/2009

/*
 element ZeroLengthContactNTS2D eleTag? -sNdNum sNode? -pNdNum pNode? -Nodes Nodes? Kn? Kt? phi?
 Description: This file contains the implementation for the ZeroLengthContactNTS2D class.
*/

#include "ZeroLengthContactNTS2D.h"
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ElementResponse.h>

static int numZeroLengthContactNTS2D = 0;
#include <elementAPI.h>

#define PI 3.141592653589793238462643383279502884197169399

void * 
OPS_ZeroLengthContactNTS2D(void) {

  if (numZeroLengthContactNTS2D == 0) {
    numZeroLengthContactNTS2D++;
    opserr << "ZeroLengthContactNTS2d - Written by Roozbeh G. Mikola and N.Sitar, UC Berkeley\n";
  }

  Element *theEle = 0;
  int numData = 0;  

  // get the ele tag
  int eleTag, sNdNum, pNdNum;
  numData = 1;

  if (OPS_GetInt(&numData, &eleTag) != 0) {
    opserr << "ZeroLengthContactNTS2D::WARNING invalied eleTag \n";
    return 0;
  }

  const char *nextString = OPS_GetString();
  if (strcmp(nextString,"-sNdNum") != 0) {
    opserr << "ZeroLengthContactNTS2D:: expecting "<<
      "- element ZeroLengthContactNTS2D eleTag? -sNdNum sNode? -pNdNum pNode? -Nodes Nodes? Kn? Kt? phi? \n" ;
    return 0;
  }

  // get the number of secondary nodes
  numData = 1;
  if (OPS_GetInt(&numData, &sNdNum) != 0) {
    opserr << "ZeroLengthContactNTS2D::WARNING invalied sNdNum \n";
    return 0;
  }

  numData = 10;
  nextString = OPS_GetString();

  if (strcmp(nextString,"-mNdNum") != 0 && strcmp(nextString,"-pNdNum") != 0) {
    opserr << "ZeroLengthContactNTS2D:: expecting "<<
      "- element ZeroLengthContactNTS2D eleTag? -sNdNum sNode? -pNdNum pNode? -Nodes Nodes? Kn? Kt? phi? \n" ;
    return 0;
  }
  
  numData = 1;
  if (OPS_GetInt(&numData, &pNdNum) != 0) {
    opserr << "ZeroLengthContactNTS2D::WARNING invalied sNdNum \n";
    return 0;
  }

  // a quick check on number of args
  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 3 + sNdNum + pNdNum) {
    opserr << "ZeroLengthContactNTS2D::WARNING too few arguments " <<
      "want - element zeroLengthContactNTS2D $tag -sNdNum $sNdNum -pNdNum $pNdNum -Nodes $Nodes $Kn $Kt $phi" ;
    return 0;
  }

  numData = 10;
  nextString = OPS_GetString();

  if (strcmp(nextString,"-Nodes") != 0) {
    opserr << "ZeroLengthContactNTS2D:: expecting "<<
      "- element ZeroLengthContactNTS2D eleTag? -sNdNum sNode? -pNdNum pNode? -Nodes Nodes? Kn? Kt? phi? \n" ;
    return 0;
  }


  // read the Nodes values
  numData = sNdNum+pNdNum;
  int *theNodeData = new int[numData];
  ID  Nodes(theNodeData, numData);

  if (OPS_GetInt(&numData, theNodeData) != 0) {
    opserr << "ZeroLengthContactNTS2D:: invalid Nodes number value for -Nodes ";
    opserr << eleTag << "- element ZeroLengthContactNTS2D eleTag? -sNdNum sNode? -pNdNum pNode? -Nodes Nodes? Kn? Kt? phi? \n" ;
    return 0;
  }

  // read the material properties
  numData = 3;
  double dData[3];
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "ZeroLengthContactNTS2D::WARNING invalid Kn,Kt or phi\n" ;
    return 0;
  }

  //
  // now we create the element and add it to the domain
  //
  
  theEle = new ZeroLengthContactNTS2D(eleTag, sNdNum, pNdNum, Nodes, dData[0], dData[1], dData[2]);
  return theEle;
}
  



//*********************************************************************
//  Full Constructor:







ZeroLengthContactNTS2D::ZeroLengthContactNTS2D(int tag, int sNdNum, int pNdNum, const ID& Nodes,
					       double Knormal, double Ktangent, double frictionAngle)
  :Element(tag,ELE_TAG_ZeroLengthContactNTS2D),
   connectedExternalNodes(sNdNum + pNdNum),
   N(6), T(6), ContactNormal(2), Ki(0), load(0)
{
  //static data
  numberNodes = sNdNum + pNdNum;
  SecondaryNodeNum = sNdNum;
  PrimaryNodeNum = pNdNum;
    
  // static data for 2D
  stiff.resize(2*numberNodes,2*numberNodes);
  resid.resize(2*numberNodes);
  zeroMatrix.resize(2*numberNodes,2*numberNodes);
  
  // allocate contact force vector
  pressure.resize(numberNodes);
  // allocate normal gap vector
  normal_gap.resize(numberNodes);
  // allocate shear gap vector
  shear_gap.resize(numberNodes);
  stored_shear_gap.resize(numberNodes);
  
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != numberNodes)
    opserr << "FATAL ZeroLength::setUp - failed to create an ID of correct size\n";
  
  // create an ID of correct size
  nodePointers = new Node* [numberNodes];
  
  // set the vectors to zero
  for(int i = 0; i < numberNodes; i++)
    {
      stored_shear_gap(i) = 0;
      shear_gap(i) = 0;
      pressure(i) = 0;
      normal_gap(i) = 0;
    }
  
  // restore node number
  for (int i = 0 ; i < numberNodes; i++) 
    connectedExternalNodes(i) = Nodes(i);
  
  // assign Kn, Kt, fc
  Kn = Knormal;
  Kt = Ktangent;
  // friction ratio fc = tan(phi)
  fc = tan (frictionAngle * PI / 180);
  
  // initialized contact flag be zero
  ContactFlag=0;
}

//null constructor
ZeroLengthContactNTS2D::ZeroLengthContactNTS2D(void)
  :Element(0,ELE_TAG_ZeroLengthContactNTS2D),
  connectedExternalNodes(numberNodes),
  N(2*numberNodes), T(2*numberNodes), Ki(0), load(0)
{
    // ensure the connectedExternalNode ID is of correct size
    if (connectedExternalNodes.Size() != numberNodes)
		opserr << "FATAL ZeroLengthContactNTS2D::ZeroLengthContactNTS2D - failed to create an ID of correct size\n";
    for (int j = 0; j < numberNodes; j++ ) 
		nodePointers[j] = 0;
}


//  Destructor:
//  delete must be invoked on any objects created by the object
//  and on the matertial object.
ZeroLengthContactNTS2D::~ZeroLengthContactNTS2D()
{
  if (load != 0) delete load;
  if (Ki != 0) delete Ki;
  //delete stored_shear_gap;
}

int
ZeroLengthContactNTS2D::getNumExternalNodes(void) const
{
    return numberNodes;
}

const ID &
ZeroLengthContactNTS2D::getExternalNodes(void)
{
  return connectedExternalNodes;
}

Node **
ZeroLengthContactNTS2D::getNodePtrs(void)
{
    return nodePointers;
}

int
ZeroLengthContactNTS2D::getNumDOF(void)
{
     return numDOF;
}

// method: setDomain()
// to set a link to the enclosing Domain and to set the node pointers.
// also determines the number of dof associated
// with the ZeroLengthContactNTS2D element
void
ZeroLengthContactNTS2D::setDomain(Domain *theDomain)
{
	// check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
		for (int j = 0; j < numberNodes; j++ )  
			nodePointers[j] = 0;
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
			opserr << "WARNING ZeroLengthContactNTS2D::setDomain() - Nd: " << Nd << " does not exist in ";
            return; 
		}
        // now determine the number of dof and the dimension
		int dofNd = nodePointers[i]->getNumberDOF();
        // if differing dof at the ends - print a warning message
        if ( dofNd != 2 ) {
			opserr << "WARNING ZeroLengthContactNTS2D::setDomain cannot handle " << dofNd << "dofs\n";
	        return;
        }
		// summation of DOFs
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
      opserr << "WARNING ZeroLengthContactNTS2D::setDomain(): Element " << this->getTag() << " has L= " << L <<
	", which is greater than the tolerance\n";
*/

}

int
ZeroLengthContactNTS2D::commitState()
{
   // need to update stick point here
    if (ContactFlag == 2 )   // slide case, update stick point
		for(int i = 0; i < numberNodes; i++) 
			stored_shear_gap(i) = shear_gap(i);

	return 0;
}

int
ZeroLengthContactNTS2D::revertToLastCommit()
{
	///////////////////////////////////////////
    // need to revert the stored_shear_gap??
    //int code=0;
    // revert state for 1d materials
    //for (int i=0; i<numMaterials1d; i++)
	//code += theMaterial1d[i]->revertToLastCommit();
    //return code;

	return 0;

}

int
ZeroLengthContactNTS2D::revertToStart()
{
	// need to rezero stored_shear_gap??
    //int code=0;
    // revert to start for 1d materials
    //for (int i=0; i<numMaterials1d; i++)
	//code += theMaterial1d[i]->revertToStart();
    //return code;
	// zero stored_shear_gap

	for(int i = 0; i < numberNodes; i++) 
		stored_shear_gap(i) = 0;
	return 0;
}

// update
// calculate stress-strain relation -- M. Frank
int
ZeroLengthContactNTS2D::update(void)
{
	return 0;
}

const Matrix &
ZeroLengthContactNTS2D::getTangentStiff(void)
{
  int tang_flag = 1 ; //get the tangent
  //zero stiffness and residual
  stiff.Zero();
  //do tangent and residual here
  formGlobalResidAndTangent( tang_flag ) ;

  return stiff ;
}

const Matrix &
ZeroLengthContactNTS2D::getInitialStiff(void)
{
  int tang_flag = 1 ; //get the tangent
  stiff.Zero();
  //do tangent and residual here
  formGlobalResidAndTangent( tang_flag ) ;
  return stiff ;
}

const Matrix &
ZeroLengthContactNTS2D::getDamp(void)
{
    // no damp
 	zeroMatrix.Zero();
	return zeroMatrix;
}

const Matrix &
ZeroLengthContactNTS2D::getMass(void)
{
    // no mass
 	zeroMatrix.Zero();
	return zeroMatrix;
}

void
ZeroLengthContactNTS2D::zeroLoad(void)
{
  // do nothing now
}

int
ZeroLengthContactNTS2D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // meaningless to addLoad to a contact !
  return 0;
}

int
ZeroLengthContactNTS2D::addInertiaLoadToUnbalance(const Vector &accel)
{
  // do nothing as element has no mass yet!
  return 0;
}

const Vector &
ZeroLengthContactNTS2D::getResistingForce()
{
  int tang_flag = 0 ; //don't get the tangent
  resid.Zero();
  formGlobalResidAndTangent( tang_flag ) ;
  //opserr<< "resid="<<resid;
  return resid ;
}

const Vector &
ZeroLengthContactNTS2D::getResistingForceIncInertia()
{
  // there is no Inertia
  int tang_flag = 0 ; //don't get the tangent
  resid.Zero();
  formGlobalResidAndTangent( tang_flag ) ;
  return  resid ;
}

int
ZeroLengthContactNTS2D::sendSelf(int commitTag, Channel &theChannel)
{
 // doing nothing here
	return 0;
}

int
ZeroLengthContactNTS2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
// doing nothing here
	return 0;
}

int
ZeroLengthContactNTS2D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{ 
 // nothing to display
    return 0;
}

void
ZeroLengthContactNTS2D::Print(OPS_Stream &s, int flag)
{
    if (flag == 0) { // print everything
		s << "Element: " << this->getTag();
	    s << " type: ZeroLengthContactNTS2D  Nodes: " << connectedExternalNodes << endln;
    } else if (flag == 1) {
		s << this->getTag() << "  ";
    }
}

Response*
ZeroLengthContactNTS2D::setResponse(const char **argv, int argc, Information &eleInformation)
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
ZeroLengthContactNTS2D::getResponse(int responseID, Information &eleInfo)
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
// determine the secondar/primary pair in contact, and setup Vectors (N,T)
int ZeroLengthContactNTS2D::contactDetect(int s, int m1, int m2, int stage)
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
	 // get current position and after trial displacement (secondary, primary1, primary2)

	 const Vector &xs = nodePointers[s]->getCrds();
	 const Vector &uxs = nodePointers[s]->getTrialDisp();
     const Vector &x1 = nodePointers[m1]->getCrds();
	 const Vector &ux1= nodePointers[m1]->getTrialDisp();
     const Vector &x2 = nodePointers[m2]->getCrds();
	 const Vector &ux2= nodePointers[m2]->getTrialDisp();
	 Vector trial_secondary = xs + uxs;
     Vector trial_primary1 = x1 + ux1;
     Vector trial_primary2 = x2 + ux2;

	 Vector diff = trial_primary2 - trial_primary1;
	 // Length of segment
     double L  = diff.Norm();
	 // tangent vector
	 Vector ContactTangent(2);
	 ContactTangent = (1/L) * (trial_primary2 - trial_primary1);
	 // normal vector
     ContactNormal(0) = -ContactTangent(1);
     ContactNormal(1) =  ContactTangent(0);
     // local coordination alpha = 0 for starting node and alpha = 1 for end node
	 double alpha = 0;
	 for (int i = 0; i < 2; i++) 
		 alpha += (1/L) * (trial_secondary(i) - trial_primary1(i)) * ContactTangent(i);
     // normal gap
	 normal_gap(s) = 0;
	 for (int i = 0; i < 2; i++) 
		 normal_gap(s) += (trial_secondary(i) - trial_primary1(i)) * ContactNormal(i);
     // Length of segment before deformation
	 diff = x2 - x1;
	 double L_bar = diff.Norm();
     // local coordination before deformation
	 double alpha_bar = 0;
	 for (int i = 0; i < 2; i++) 
		 alpha_bar += (1/L_bar) * (xs(i) - x1(i)) * ContactTangent(i);
     // shear deformation given before and after deformation along segment
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

void  ZeroLengthContactNTS2D::formLocalResidAndTangent( int tang_flag , int secondary, int primary1, int primary2, int stage)
{
	// trial frictional force vectors (in local coordinate)
    double t_trial;
    double TtrNorm;
    // Coulomb friction law surface
	double Phi;
	int i, j;

	for (i = 0; i < numberNodes; i++) 
		pressure(i) = 0;

    t_trial=0;

	//int IsContact;
	// detect contact and set flag
    ContactFlag = contactDetect(secondary,primary1,primary2, stage);

	if (ContactFlag == 1) // contacted
	{
		// create a vector for converting local matrix to global
        int loctoglob[6];
        loctoglob[0] = (2 * secondary); loctoglob[1] = (2 * secondary) + 1; 
	    loctoglob[2] = (2 * primary1); loctoglob[3] = (2 * primary1) + 1; 
	    loctoglob[4] = (2 * primary2); loctoglob[5] = (2 * primary2) + 1; 

		// contact presure;
	    pressure(secondary) = Kn * normal_gap(secondary);  // pressure is positive if in contact

		double ng = normal_gap(secondary);

		t_trial = Kt * (shear_gap(secondary) - stored_shear_gap(secondary));  // trial shear force

        // Coulomb friction law, trial state
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

void  ZeroLengthContactNTS2D::formGlobalResidAndTangent( int tang_flag )
{
	// in the first loop the node to node contact will not be considered 
	// in order to prevent node-node contact duplication 
	// but on contrary in the second loop the node to node contact 
	// will be considered and this can be controlled by "stage = 0 or 1"

	// loop over sedondary nodes and find the nodes 
    // which are in contact with primary's segments
	for (int i = 0 ; i < SecondaryNodeNum; i++) {
		for (int j = SecondaryNodeNum ; j < SecondaryNodeNum + PrimaryNodeNum - 1; j++) {
			formLocalResidAndTangent( tang_flag, i, j, j+1 , 0);  // stage = 0 //
        } // endfor j
    } // endfor i

    // loop over primary nodes and find the nodes 
    // which are in contact with secondary's segments
	for (int i = SecondaryNodeNum ; i < SecondaryNodeNum + PrimaryNodeNum; i++) {
		for (int j = 0 ; j < SecondaryNodeNum - 1; j++) {
            formLocalResidAndTangent( tang_flag, i, j, j+1 , 1);  // stage = 1 //
        } // endfor j
    } // endfor i
}
