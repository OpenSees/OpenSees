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

// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLengthContact2D.cpp,v $
// $Revision: 1.4 $
// $Date: 2008-06-13 21:01:03 $

// Written: Gang Wang  (wang@ce.berkeley.edu)
//          Prof. Nicholas Sitar (nsitar@ce.berkeley.edu)
//
// Created: 27/08/2003

/*
 element zeroLengthContact2D  $eleID $sNdID $mNdID $Kn $Kt $fs -normal $Nx $Ny
 Description: This file contains the implementation for the ZeroLengthContact2D class.
*/

#include "ZeroLengthContact2D.h"
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
#include <elementAPI.h>

//static data
const int ZeroLengthContact2D::numberNodes = 2 ;

// static data for 2D
Matrix  ZeroLengthContact2D::stiff4(2*numberNodes,2*numberNodes) ;
Vector  ZeroLengthContact2D::resid4(2*numberNodes) ;

Matrix  ZeroLengthContact2D::stiff6(3*numberNodes,3*numberNodes) ;
Vector  ZeroLengthContact2D::resid6(3*numberNodes) ;

void* OPS_ZeroLengthContact2D()
{
    if (OPS_GetNumRemainingInputArgs() < 9) {
	opserr << "ZeroLengthContact2D::WARNING too few arguments " <<
	    "want - element ZeroLengthContact2D eleTag? iNode? jNode? Kn? Kt? fs? -normal Nx? Ny?" ;
	return 0;
    }

    // eleTag, iNode, jNode;
    int idata[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata,idata) < 0) {
	opserr<<"WARNING: invalid integer inputs\n";
	return 0;
    }

    // Kn, Kt, fs
    double data[3];
    numdata = 3;
    if (OPS_GetDoubleInput(&numdata,data) < 0) {
	opserr<<"WARNING: invalid double inputs\n";
	return 0;
    }

    const char* type = OPS_GetString();
    if (strcmp(type,"-normal") != 0) {
	opserr << "ZeroLengthContact2D:: expecting "<<
	    "- element ZeroLengthContact2D eleTag? iNode? jNode? Kn? Kt? fs? -normal Nx? Ny? \n" ;
	return 0;
    }

    Vector normaldir(2);
    numdata = 2;
    if (OPS_GetDoubleInput(&numdata, &normaldir(0)) < 0) {
	opserr<<"WARNING: invalid double inputs\n";
	return 0;
    }

    return new ZeroLengthContact2D(idata[0],idata[1],idata[2],data[0],data[1],data[2],normaldir);

}


//*********************************************************************
//  Full Constructor:

ZeroLengthContact2D::ZeroLengthContact2D(int tag,
					 int Nd1, int Nd2,
					 double Knormal, double Ktangent,
					 double frictionRatio,  const Vector& normal )
  :Element(tag,ELE_TAG_ZeroLengthContact2D),
   connectedExternalNodes(numberNodes),
   N(2*numberNodes), T(2*numberNodes), ContactNormal(2)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)
      opserr << "FATAL ZeroLength::setUp - failed to create an ID of correct size\n";
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;

	// assign Kn, Kt, fs
	Kn = Knormal;
	Kt = Ktangent;
	fs = frictionRatio;

	// assign outward contact normal of primary block
 	ContactNormal(0) = normal(0)/normal.Norm();
	ContactNormal(1) = normal(1)/normal.Norm();

    // set stick point cords in LOCAL basis
	stickPt = 0;

	// initialized contact flag be zero
	ContactFlag=0;

	gap_n  = 0 ;
	lambda =0;   // add for augmented lagrange
	pressure=0;  // add for augmented lagrange
}



//null constructor

ZeroLengthContact2D::ZeroLengthContact2D(void)
  :Element(0,ELE_TAG_ZeroLengthContact2D),
  connectedExternalNodes(numberNodes),
  N(2*numberNodes), T(2*numberNodes), ContactNormal(2)
{

  //opserr<<this->getTag()<< " new ZeroLengthContact2D::null constructor" <<endln;

    // ensure the connectedExternalNode ID is of correct size
    if (connectedExternalNodes.Size() != 2)
      opserr << "FATAL ZeroLengthContact2D::ZeroLengthContact2D - failed to create an ID of correct size\n";
    for (int j=0; j<numberNodes; j++ ) {
      nodePointers[j] = 0;
	}
}


//  Destructor:
//  delete must be invoked on any objects created by the object
//  and on the matertial object.
ZeroLengthContact2D::~ZeroLengthContact2D()
{

  //opserr<<this->getTag()<<" ZeroLengthContact2D::destructor" <<endln;


}


int
ZeroLengthContact2D::getNumExternalNodes(void) const
{

	//opserr<<this->getTag()<<" ZeroLengthContact2D::getNumExternalNodes" <<endln;
    return 2;
}


const ID &
ZeroLengthContact2D::getExternalNodes(void)
{
  //opserr<<this->getTag()<<" ZeroLengthContact2D::getExternalNodes" <<endln;
  return connectedExternalNodes;
}



Node **
ZeroLengthContact2D::getNodePtrs(void)
{
    return nodePointers;
}

int
ZeroLengthContact2D::getNumDOF(void)
{
     return numDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the ZeroLengthContact2D element


void
ZeroLengthContact2D::setDomain(Domain *theDomain)
{

    //opserr<<this->getTag()<< " ZeroLengthContact2D::setDomain" <<endln;


    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	nodePointers[0] = 0;
	nodePointers[1] = 0;
	return;
    }

    // set default values for error conditions

    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    nodePointers[0] = theDomain->getNode(Nd1);
    nodePointers[1] = theDomain->getNode(Nd2);

    // if can't find both - send a warning message
    if ( nodePointers[0] == 0 || nodePointers[1] == 0 ) {
      if (nodePointers[0] == 0)
        opserr << "WARNING ZeroLengthContact2D::setDomain() - Nd1: " << Nd1 << " does not exist in ";
      else
        opserr << "WARNING ZeroLengthContact2D::setDomain() - Nd2: " << Nd2 << " does not exist in ";

      //opserr << "model for ZeroLengthContact2D ele: " << this->getTag() << endln;

      return;
    }

    // now determine the number of dof and the dimension
    int dofNd1 = nodePointers[0]->getNumberDOF();
    int dofNd2 = nodePointers[1]->getNumberDOF();

    // if differing dof at the ends - print a warning message
    if ( dofNd1 != dofNd2 ) {
      opserr << "WARNING ZeroLengthContact2D::setDomain(): nodes " << Nd1 << " and " << Nd2 <<
	"have differing dof at ends for ZeroLengthContact2D " << this->getTag() << endln;
      return;
    }

    // Check that length is zero within tolerance
    const Vector &end1Crd = nodePointers[0]->getCrds();
    const Vector &end2Crd = nodePointers[1]->getCrds();
    Vector diff = end1Crd - end2Crd;
    double L  = diff.Norm();
    double v1 = end1Crd.Norm();
    double v2 = end2Crd.Norm();
    double vm;

    vm = (v1<v2) ? v2 : v1;


    if (L > LENTOL*vm)
      opserr << "WARNING ZeroLengthContact2D::setDomain(): Element " << this->getTag() << " has L= " << L <<
	", which is greater than the tolerance\n";

    // call the base class method
    this->DomainComponent::setDomain(theDomain);

	if (dofNd1 == 2 && dofNd2 == 2) {
	numDOF = 4;
	}
	else if (dofNd1 == 3 && dofNd2 == 3) {
	  stiff6.Zero();
	  resid6.Zero();
	numDOF = 6;
	}	
    else {
    opserr << "WARNING ZeroLengthContact2D::setDomain cannot handle " << dofNd1 <<
 	"dofs at nodes in " << dofNd1 << " d problem\n";
     return;
    }

}

int
ZeroLengthContact2D::commitState()

{

   // need to update stick point here



    if (ContactFlag == 2 )   // slide case, update stick point

    stickPt=xi;



    gap_n = gap;



	////////////////////////////////////

	// initialize lagrange multiplier zero for next iterations.

	// lambda = pressure;          // using for augmented lagrange

	lambda = 0;                 // using penalty method only

	////////////////////////////////////



	return 0;

}



int

ZeroLengthContact2D::revertToLastCommit()

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

ZeroLengthContact2D::revertToStart()

{



	// need to rezero stickPoint??





    //int code=0;

    // revert to start for 1d materials

    //for (int i=0; i<numMaterials1d; i++)

	//code += theMaterial1d[i]->revertToStart();

    //return code;





	// zero stickPt

	stickPt=0;

	return 0;

}





// update

// calculate stress-strain relation -- M. Frank



int

ZeroLengthContact2D::update(void)

{

	//////////////////////////////////////////////////////////

    //opserr<<"contact2D update " << endln;

	//lambda   = 990 ;  // add for augmented lagrange

	//opserr<<"contact2D lambda = " << lambda <<endln;

	///////////////////////////////////////////////////////////

	return 0;

}





const Matrix &

ZeroLengthContact2D::getTangentStiff(void)

{



  //opserr<< this->getTag()<< " ZeroLengthContact2D::getTangentStiff()" <<endln;



  int tang_flag = 1 ; //get the tangent



  //do tangent and residual here

  formResidAndTangent( tang_flag ) ;



  //opserr<< stiff ;


  if (numDOF == 4)
    return stiff4;
  else {
    stiff6(0,0) = stiff4(0,0);
    stiff6(1,0) = stiff4(1,0);
    stiff6(0,1) = stiff4(0,1);
    stiff6(1,1) = stiff4(1,1);

    stiff6(3,3) = stiff4(2,2);
    stiff6(4,3) = stiff4(3,2);
    stiff6(3,4) = stiff4(2,3);
    stiff6(4,4) = stiff4(3,3);    

    stiff6(0,3) = stiff4(0,2);
    stiff6(1,3) = stiff4(1,2);
    stiff6(0,4) = stiff4(0,3);
    stiff6(1,4) = stiff4(1,3);

    stiff6(3,0) = stiff4(2,0);
    stiff6(4,0) = stiff4(3,0);
    stiff6(3,1) = stiff4(2,1);
    stiff6(4,1) = stiff4(3,1);

    return stiff6;
  }



}





const Matrix &

ZeroLengthContact2D::getInitialStiff(void)

{



  //opserr<<this->getTag()<< " ZeroLengthContact2D::getInitialStiff()" <<endln;



  int tang_flag = 1 ; //get the tangent



  //do tangent and residual here

  formResidAndTangent( tang_flag ) ;


  if (numDOF == 4)
    return stiff4;
  else {
    stiff6(0,0) = stiff4(0,0);
    stiff6(1,0) = stiff4(1,0);
    stiff6(0,1) = stiff4(0,1);
    stiff6(1,1) = stiff4(1,1);

    stiff6(3,3) = stiff4(2,2);
    stiff6(4,3) = stiff4(3,2);
    stiff6(3,4) = stiff4(2,3);
    stiff6(4,4) = stiff4(3,3);    

    stiff6(0,3) = stiff4(0,2);
    stiff6(1,3) = stiff4(1,2);
    stiff6(0,4) = stiff4(0,3);
    stiff6(1,4) = stiff4(1,3);

    stiff6(3,0) = stiff4(2,0);
    stiff6(4,0) = stiff4(3,0);
    stiff6(3,1) = stiff4(2,1);
    stiff6(4,1) = stiff4(3,1);

    return stiff6;
  }

}





const Matrix &

ZeroLengthContact2D::getDamp(void)

{

    // no damp

  if (numDOF == 4) {
    stiff4.Zero();
    return stiff4;
  } else {
    stiff6.Zero();
    return stiff6;
  }
}





const Matrix &

ZeroLengthContact2D::getMass(void)

{

    // no mass
  if (numDOF == 4) {
    stiff4.Zero();
    return stiff4;
  } else {
    stiff6.Zero();
    return stiff6;
  }
}




void
ZeroLengthContact2D::zeroLoad(void)
{
  // does nothing now

}

int
ZeroLengthContact2D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // meaningless to addLoad to a contact !
  return 0;
}

int
ZeroLengthContact2D::addInertiaLoadToUnbalance(const Vector &accel)
{
  // does nothing as element has no mass yet!
  return 0;
}


const Vector &
ZeroLengthContact2D::getResistingForce()
{

  //opserr<< this->getTag()<< " ZeroLengthContact2D::getResistingForce()" <<endln;

  int tang_flag = 0 ; //don't get the tangent
  formResidAndTangent( tang_flag ) ;

  //opserr<< "resid="<<resid;

  if (numDOF == 4)
    return resid4;
  else {
    resid6(0) = resid4(0);
    resid6(1) = resid4(1);
    resid6(3) = resid4(2);
    resid6(4) = resid4(3);
    return resid6;
  }
}


const Vector &
ZeroLengthContact2D::getResistingForceIncInertia()
{
    // there is no Inertia
  //opserr<< this->getTag()<< " ZeroLengthContact2D::getResistingForceIncInertia()" <<endln;

  int tang_flag = 0 ; //don't get the tangent
  formResidAndTangent( tang_flag ) ;

  //opserr<< resid;

  if (numDOF == 4)
    return resid4;
  else {
    resid6(0) = resid4(0);
    resid6(1) = resid4(1);
    resid6(3) = resid4(2);
    resid6(4) = resid4(3);
    return resid6;
  }
}


int
ZeroLengthContact2D::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  int dbTag = this->getDbTag();

  static ID idData(5);
  idData(0) = this->getTag();
  idData(1) = connectedExternalNodes(0);
  idData(2) = connectedExternalNodes(1);  
  idData(3) = numDOF;
  idData(4) = ContactFlag;

  res += theChannel.sendID(dbTag, commitTag, idData);
  if (res < 0) {
    opserr << "ZeroLengthContact2D::sendSelf -- failed to send ID data" << endln;
    return res;
  }

  static Vector data(10);
  data(0) = pressure;
  data(1) = lambda;
  data(2) = t1;
  data(3) = t2;
  data(4) = gap_n;
  data(5) = Kn;
  data(6) = Kt;
  data(7) = fs;
  data(8) = stickPt;
  data(9) = xi;
  
  res += theChannel.sendVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "ZeroLengthContact2D::sendSelf -- failed to send Vector data" << endln;
    return res;
  }
  
  return res;
}

int
ZeroLengthContact2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dbTag = this->getDbTag();

  static ID idData(5);
  res += theChannel.recvID(dbTag, commitTag, idData);
  if (res < 0) {
    opserr << "ZeroLengthContact2D::recvSelf -- failed to receive ID data" << endln;
    return res;
  }
  
  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);
  numDOF = idData(3);
  ContactFlag = idData(4);

  static Vector data(10);
  res += theChannel.recvVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "ZeroLengthContact2D::recvSelf -- failed to receive Vector data" << endln;
    return res;
  }

  pressure = data(0);
  lambda = data(1);
  t1 = data(2);
  t2 = data(3);
  gap_n = data(4);
  gap = gap_n;
  Kn = data(5);
  Kt = data(6);
  fs = data(7);
  stickPt = data(8);
  xi = data(9);
  
  return res;
 
}


int
ZeroLengthContact2D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{ // nothing to display
    return 0;
}


void
ZeroLengthContact2D::Print(OPS_Stream &s, int flag)
{
    if (flag == 0) { // print everything
	s << "Element: " << this->getTag();
	s << " type: ZeroLengthContact2D  iNode: " << connectedExternalNodes(0);
	s << " jNode: " << connectedExternalNodes(1) << endln;
    } else if (flag == 1) {
	s << this->getTag() << "  ";
    }

}

Response*
ZeroLengthContact2D::setResponse(const char **argv, int argc, OPS_Stream &output)
{
     if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
       return new ElementResponse(this, 1, Vector(numDOF));

     // tangent stiffness matrix
     else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0)
       return new ElementResponse(this, 2, Matrix(numDOF,numDOF));

	 // contact pressure
     else if (strcmp(argv[0],"pressure")== 0)
	{//opserr<<"Contact2DsetResponse p="<<this->pressure <<endln;
     return new ElementResponse(this, 3, pressure);
	 }

     else if (strcmp(argv[0],"gap")== 0)
     return new ElementResponse(this, 4, gap);

  	else
	  return Element::setResponse(argv, argc, output);

}


int
ZeroLengthContact2D::getResponse(int responseID, Information &eleInfo)
{
 if (responseID == 1)
	 return eleInfo.setVector(this->getResistingForce());
 else if (responseID == 2)
	 return eleInfo.setMatrix(this->getTangentStiff());
 else if (responseID == 3)
 {//opserr<<"Contact2D getResponse p="<<this->pressure<<endln;
 return eleInfo.setDouble(this->pressure);
 }
 else if (responseID == 4)
  return eleInfo.setDouble(this->gap);

 else
   return Element::getResponse(responseID, eleInfo);
}


// Private methods

// determine the secondary/primary pair in contact, and setup Vectors (N,T1,T2)
 int ZeroLengthContact2D::contactDetect(void)
 {
  	//opserr<< this->getTag()<< " ZeroLengthContact2D::contactDetect" <<endln;


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
	  // gap = (U_primary-U_secondary) \dot (ContactNormal),
	  // defines overlapped normal distance, always keep positive (+) when contacted
      ///
      // get current trial position
   //const Vector &U_secondary = nodePointers[0]->getCrds() + nodePointers[0]->getTrialDisp();
   //const Vector &U_primary= nodePointers[1]->getCrds() + nodePointers[1]->getTrialDisp();
   Vector U_secondary(2);
   Vector U_primary(2);   
   const Vector &XI = nodePointers[0]->getCrds();
   const Vector &XJ = nodePointers[1]->getCrds();
   const Vector &UI = nodePointers[0]->getTrialDisp();
   const Vector &UJ = nodePointers[1]->getTrialDisp();
   U_secondary(0) = XI(0) + UI(0);
   U_secondary(1) = XI(1) + UI(1);
   U_primary(0) = XJ(0) + UJ(0);
   U_primary(1) = XJ(1) + UJ(1);
   
	       gap=0;
		   int i;
		   for (i=0; i<2; i++){
				   gap += (U_primary(i)-U_secondary(i))* ContactNormal(i);
		   }
       //*////////////////////////////// for transient gap ///////////////////////////////

	  // we have another way to define the gap, can replace previous code block if want
      /*/////////////// for dynamic gap //////////////////////
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

      ///////////////// for dynamic gap //////////////////////*/




			if (gap < 0)   // Not in contact
			   return 0;
			else
			{
				N(0)   =  ContactNormal(0);
				N(1)   =  ContactNormal(1);
				N(2)   = -N(0) ;
				N(3)   = -N(1);

				T(0)  =    N(1);
				T(1)  =   -N(0);
				T(2)  =   -T(0);
				T(3)  =   -T(1);

 			   return 1;
			}
  }


void  ZeroLengthContact2D::formResidAndTangent( int tang_flag )
{

	//opserr<<this->getTag()<< " ZeroLengthContact2D:: formResidAndTangent()" <<endln;


	// trial displacement vectors
 	static Vector DispTrialS(2); // trial disp for secondary node
	static Vector DispTrialP(2); // trial disp for primary node
	// trial frictional force vectors (in local coordinate)
    double t_trial;
    double TtrNorm;

    // Coulomb friction law surface
	double Phi;

    int i, j;

    //zero stiffness and residual
    Matrix &stiff = stiff4;
    Vector &resid = resid4;
    stiff.Zero( ) ;
    resid.Zero( ) ;

	pressure = 0;
    t_trial=0;

	//int IsContact;
	// detect contact and set flag
    ContactFlag = contactDetect();

	//opserr<<this->getTag()<< " ZeroLengthContact2D::ContactFlag=" << ContactFlag<<endln;

	if (ContactFlag == 1) // contacted
	//if (gap >= 0) // contacted
    //if ((lambda + Kn*gap) >= 0)  // changed for augmented lagrange
	{
       // contact presure;
	    pressure = Kn*gap ;  // pressure is positive if in contact
       // pressure = Kn*gap + lambda;  // changed for augmented lagrange

	    //DispTrialS=nodePointers[0]->getTrialDisp();
	    //DispTrialP=nodePointers[1]->getTrialDisp();
	    const Vector &UI = nodePointers[0]->getTrialDisp();
	    const Vector &UJ = nodePointers[1]->getTrialDisp();	    
	    DispTrialS(0) = UI(0);
	    DispTrialS(1) = UI(1);
	    DispTrialP(0) = UJ(0);
	    DispTrialP(1) = UJ(1);	    
	    
        //opserr<<"DispTrialS " << DispTrialS;
        //opserr<<"DispTrialP " << DispTrialP;

       //nodal displacements
        double ul[4];

		ul[0]=DispTrialS(0);
		ul[1]=DispTrialS(1);
 		ul[2]=DispTrialP(0);
		ul[3]=DispTrialP(1);

		t_trial = 0;
        xi=0;


		// relative slide displacement
		// xi = T_tran * u    eq. (3.5)
		for (i=0; i<4; i++){
			xi  += T (i)*ul[i];
		}

		//for (i=0; i<2; i++){ t_trial(i)=Kt * (xi(i)-stickPt(i));}  //3D

		 t_trial=Kt*(xi-stickPt);  // trial shear force
         /// t_trial=Kt*(xi);  // no update of stickPt, updated Jan 26, 2004

        // Coulomb friction law, trial state

		//TtrNorm=t_trial.Norm();

		TtrNorm=sqrt(t_trial*t_trial);

		Phi = TtrNorm - fs * pressure;


		if (Phi <= 0 ) { // stick case
  		//opserr<< "stick ...." << endln;

 			if ( tang_flag == 1 ) {
		    // stiff
				for (i=0; i<4; i++) {
					for (j=0; j<4; j++) {
						//stiff(i,j) = Kn*(N(i)*N(j)) + Kt*(T(i)*T1(j)+T2(i)*T2(j));// 3D
                        stiff(i,j) =   Kn*(N(i)*N(j)) + Kt*(T(i)*T(j));	//2D
					}
				}
			} //endif tang_flag
			// force
			for (i=0; i<4; i++)
				 resid(i)= (-1*pressure)*N(i) + t_trial*T(i);    //2D
			//	 resid(i)= (-1*pressure)*N(i) + t_trial(0)*T1(i) + t_trial(1)*T2(i) ;%3D

		} // end if stick
		else {           // slide case, non-symmetric stiff
            ContactFlag=2;  // set the contactFlag for sliding

		 	//opserr<< "sliding ...." << endln;


            if ( tang_flag == 1 ) {
				// stiff
				for (i=0; i<4; i++) {
					for (j=0; j<4; j++) {
					// 3D
					//	define: Pt1=t_trial(0)/TtrNorm;
				    //  define: Pt2=t_trial(1)/TtrNorm;
					//	define: C1=fs*Kn;
                    //  C2 term will be zero in two dimensional formulation
					// stiff(i,j) = Kn*(N(i)*N(j)) - C1*(Pt1*T1(i)*N(j)+Pt2*T2(i)*N(j))
					//	    + C2*( (1-Pt1*Pt1)*T1(i)*T1(j) -    Pt1*Pt2 *T1(i)*T2(j)
					//      - Pt1*Pt2 *T2(i)*T1(j) + (1-Pt1*Pt2)*T2(i)*T2(j)  );  //3D

					// 2D                 ???? - or + ????
					stiff(i,j) = Kn*(N(i)*N(j)) - fs*Kn* (t_trial/TtrNorm)*T(i)*N(j); //2D
					} //endfor i
				} //endfor j
			} // endif tang_flag

			// force
			double shear=fs*pressure* (t_trial/TtrNorm);

			for (i=0; i<4; i++) {
				resid(i) = (-1*pressure)*N(i) + shear *T (i) ;      //2D
			//	resid(i) = (-1*pressure)*N(i) + t1*T1(i)+t2*T2(i) ; //3D
			}
		} //endif slide

	}  // endif ContactFlag==1



	//opserr<<"gap=" << gap <<endln;
    //opserr<<"pressure= "<<pressure <<endln;
    //opserr<<"lambda=   "<<lambda <<endln;
	//opserr<<"t_trial= "<<t_trial <<endln;
    //opserr<<"stickPt= "<<stickPt <<endln;
	//opserr<<"residue= "<<resid <<endln;

	// for NOT contact, do nothing, stiff and resid are zeroes




    /* my notes:
	   the direction of residual force is always confusing ...
	      R=KU-Fext
	   is defined as the resisting force that could provided by the element


       Let p,shear be absolute value of pressure and shear force, always positive

       thus,

                                Rx(1)=p*(-n)
                                 ||
							    \||/
						      ___\/____
						     /         \       /\
             Rx(1)     ---\ /    (1)    \     /||\n    Note: (t,n) follows RightHand rule
	         =shear*t  ---/ \   secondary   /      ||
							 \_________/       ||_____\t
						-----------------------*------/
                        |                      |
                        |                      |
						|    (2) Primary        |/---- Rx(2) = shear*(-t)
						|                      |\----
						------------------------
                                  /\
								 /||\
                                  ||
								 Rx(2)=pn

     Denote :  N={n; -n}; T={t; -t}; arrange resid={R(1); R(2)}

     Finally,  resid(i) = (- p)*N(i) + shear *T (i) ;


  */

}


