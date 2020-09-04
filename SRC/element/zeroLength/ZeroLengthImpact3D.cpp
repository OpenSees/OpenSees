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
                                                                        
// Developed by: Prof. Arash E. Zaghi & Majid Cashany, University of Connecticut (UConn) Copyright 2012
//   Based on ZereLengthContact element by Gang Wang.
//
// What: "@(#) ZeroLengthImpact3D.C, revA"


// we specify what header files we need
#include "ZeroLengthImpact3D.h"
#include <elementAPI.h>
#include <G3Globals.h>

#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// initialise the class wide variables

static int numMyZeroLengthImpact3D = 0;

void *
OPS_ZeroLengthImpact3D(void)
{
  // print out a message about who wrote this element & any copyright info wanted
  if (numMyZeroLengthImpact3D == 0) {
    opserr << "Using ZeroLengthImpact3D element - Developed by Prof. Arash E. Zaghi & Majid Cashany, University of Connecticut (UConn) Copyright 2012 - Use at your Own Peril\n";
    numMyZeroLengthImpact3D++;
  }

  // get the id and end nodes 
  int iData[4];
  double dData[7];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[0]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D tag\n";
    return 0;
  }

  int eleTag = iData[0];

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D 1st node " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[2]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D 2nd node " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[3]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D direction " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[0]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D initial gap input " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D frictionRatio " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[2]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D Ktangent " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[3]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D Knormal " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[4]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D Kn2 Input " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[5]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D Delta_y Input " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[6]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D cohesion " << eleTag << endln;
    return 0;
  }

  /*
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[7]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D origin X " << eleTag << endln;
    return 0;
  }
  
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[8]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D origin Y " << eleTag << endln;
    return 0;
  }
  */

  // now create the ZeroLengthImpact3D and add it to the Domain

  Element *theZeroLengthImpact3D = new ZeroLengthImpact3D(eleTag, iData[1], iData[2], iData[3], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6]);

  if (theZeroLengthImpact3D == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    //delete theMaterial;
    return 0;
  }

  return theZeroLengthImpact3D;
}

//static data
const int ZeroLengthImpact3D::numberNodes = 2 ;
// static data for 3D
Matrix  ZeroLengthImpact3D::stiff(3*numberNodes,3*numberNodes) ;
Vector  ZeroLengthImpact3D::resid(3*numberNodes) ;
Matrix  ZeroLengthImpact3D::zeroMatrix(3*numberNodes,3*numberNodes) ;

// typical constructor
ZeroLengthImpact3D::ZeroLengthImpact3D(int tag, 
		   int Nd1, int Nd2, 
		   int direction, 
		   double initGapInput, double frictionRatio, double Ktangent, 
		   double Knormal, double Kn2Input, double Delta_yInput,
		   double c )
:Element(tag, ELE_TAG_ZeroLengthContact3D),     
 connectedExternalNodes(numberNodes),
 directionID(direction), N(3*numberNodes), T1(3*numberNodes), T2(3*numberNodes),
 Ki(0), load(0), origin(2), stickPt(2), xi(2)
{	
  if ( direction < 0 || direction > 3 ) {
    opserr << "WARNING ZeroLengthContact3D:incorrect direction, set to 0\n";
    directionID = 0;
  }
  
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)
    opserr << "FATAL ZeroLength::setUp - failed to create an ID of correct size\n";    
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  // assign Kn, Kt, fs, cohesion
  Kn = Knormal;
  Kt = Ktangent;
  fs = frictionRatio;
  cohesion = c; 
  
  //set origin coord for circular contact, (0,0) by default
  //////////////////////////////////////////////////// removing "origX" and "origY" from input arguments: 
  //origin(0) = origX;
  origin(0) = 0.0;
  //origin(1) = origY;
  origin(1) = 0.0;
  ////////////////////////////////////////////////////
  // set stick point cords in LOCAL basis
  stickPt(0)= 0;
  stickPt(1)= 0;
  
  // initialized contact flag be zero
  ContactFlag=0;
  
  gap_n = 0 ; 

  initGap = initGapInput ;
  Kn1 = Knormal;
  Kn2 = Kn2Input; 
  Delta_y = Delta_yInput;
}

// constructor which should be invoked by an FE_ObjectBroker only
ZeroLengthImpact3D::ZeroLengthImpact3D()
  :Element(0, ELE_TAG_ZeroLengthContact3D),
   connectedExternalNodes(numberNodes),
   N(3*numberNodes), T1(3*numberNodes), T2(3*numberNodes),
   Ki(0), load(0), origin(2), stickPt(2),  xi(2)
{
  // ensure the connectedExternalNode ID is of correct size 
  if (connectedExternalNodes.Size() != 2)
    opserr << "FATAL ZeroLengthContact3D::ZeroLengthContact3D - failed to create an ID of correct size\n";
  for (int j=0; j<numberNodes; j++ ) {
    nodePointers[j] = 0;
  }
}

//  destructor - provided to clean up any memory
ZeroLengthImpact3D::~ZeroLengthImpact3D()
{
    // clean up the memory associated with the element, this is
    // memory the ZeroLengthImpact3D objects allocates and memory allocated 
    // by other objects that the ZeroLengthImpact3D object is responsible for 
    // cleaning up, i.e. the MaterialObject.
    
    if (load != 0)
    delete load;
  
    if (Ki != 0)
    delete Ki;
}

int
ZeroLengthImpact3D::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
ZeroLengthImpact3D::getExternalNodes(void) 
{
  return connectedExternalNodes;
}

Node **
ZeroLengthImpact3D::getNodePtrs(void) 
{
  return nodePointers;
}

int
ZeroLengthImpact3D::getNumDOF(void) {
    return numDOF;
}

// method: setDomain()
//    to set a link to the enclosing Domain, ensure nodes exist in Domain
//    and set pointers to these nodes, also determines the length and 
//    transformation Matrix.
void
ZeroLengthImpact3D::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	nodePointers[0] = 0;
	nodePointers[1] = 0;
	return;
    }
    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    nodePointers[0] = theDomain->getNode(Nd1);
    nodePointers[1] = theDomain->getNode(Nd2);	

    // if can't find both - send a warning message
    if ( nodePointers[0] == 0 || nodePointers[1] == 0 ) {
      if (nodePointers[0] == 0) 
        opserr << "WARNING ZeroLengthContact3D::setDomain() - Nd1: " << Nd1 << " does not exist in ";
      else
        opserr << "WARNING ZeroLengthContact3D::setDomain() - Nd2: " << Nd2 << " does not exist in ";
      return;
    }

    // now determine the number of dof and the dimension    
    int dofNd1 = nodePointers[0]->getNumberDOF();
    int dofNd2 = nodePointers[1]->getNumberDOF();	

    // if differing dof at the ends - print a warning message
    if ( dofNd1 != dofNd2 ) {
      opserr << "WARNING ZeroLengthContact3D::setDomain(): nodes " << Nd1 << " and " << Nd2 <<
	"have differing dof at ends for ZeroLengthContact3D " << this->getTag() << endln;
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
      opserr << "WARNING ZeroLengthContact3D::setDomain(): Element " << this->getTag() << " has L= " << L << 
	", which is greater than the tolerance\n";
        
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
	if (dofNd1 == 3 && dofNd2 == 3) {
	numDOF = 6;	
	}
    else {
    opserr << "WARNING ZeroLengthContact3D::setDomain cannot handle " << dofNd1 << 
 	"dofs at nodes, can only handle 3\n"; 
     return;
    }
}   	 


int
ZeroLengthImpact3D::commitState()
{
    // need to update stick point here 
  if (ContactFlag == 2 )   // update stick point only for slide case
    stickPt=xi; 

  // update gap for "dynamic gap" method
    gap_n = gap; 


	pressC = pressT;
	gapC = gapT;



	return 0;
}

int
ZeroLengthImpact3D::revertToLastCommit()
{
    ///////////////////////////////////////////
    // need to revert the stickPoint??
	xi=stickPt;


	
	pressT = pressC  ;
	gapT = gapC  ;



	return 0;
}

int
ZeroLengthImpact3D::revertToStart()
{
    // need to rezero stickPoint??
	stickPt.Zero();  


	
    gapC = 0.0;
	pressC = 0.0;
    gapT = 0.0;
	pressT = 0.0;
	tangentT = Kn1;



	return 0;
}

//int
//ZeroLengthImpact3D::update()
//{
//  // determine the current strain given trial displacements at nodes
//  double strain = this->computeCurrentStrain();
//
//  // set the strain in the materials
//  theMaterial->setTrialStrain(strain);
//
//  return 0;
//}


const Matrix &
ZeroLengthImpact3D::getTangentStiff(void)
{
    int tang_flag = 1 ; //get the tangent 

    //do tangent and residual here
    formResidAndTangent( tang_flag ) ;  

    return stiff ;
}

const Matrix &
ZeroLengthImpact3D::getInitialStiff(void)
{
    int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}

const Matrix &
ZeroLengthImpact3D::getDamp(void)
{
    // no mass 
 	zeroMatrix.Zero(); 
	return zeroMatrix;
}

const Matrix &
ZeroLengthImpact3D::getMass(void)
{
    // no mass 
 	zeroMatrix.Zero(); 
	return zeroMatrix;
}

void 
ZeroLengthImpact3D::zeroLoad(void)
{
  // does nothing now
}

int 
ZeroLengthImpact3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // meaningless to addLoad to a contact !
  return 0;
}

int 
ZeroLengthImpact3D::addInertiaLoadToUnbalance(const Vector &accel)
{
    // does nothing as element has no mass yet!
    return 0;
}

const Vector &
ZeroLengthImpact3D::getResistingForce()
{	
    int tang_flag = 0 ; //don't get the tangent
    formResidAndTangent( tang_flag ) ;
    
    return resid ;
}

const Vector &
ZeroLengthImpact3D::getResistingForceIncInertia()
{	
    // there is no Inertia 
    int tang_flag = 0 ; //don't get the tangent
    formResidAndTangent( tang_flag ) ;
   
   return  resid ;   
}

int
ZeroLengthImpact3D::sendSelf(int commitTag, Channel &theChannel)
{
    // doing nothing here
    return 0;
}

int
ZeroLengthImpact3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // doing nothing here
    return 0;
}

int
ZeroLengthImpact3D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // nothing to display
    return 0;
}

void
ZeroLengthImpact3D::Print(OPS_Stream &s, int flag)
{
    if (flag == 0) { // print everything
	s << "Element: " << this->getTag(); 
	s << " type: ZeroLengthContact3D  iNode: " << connectedExternalNodes(0);
	s << " jNode: " << connectedExternalNodes(1) << endln;
    } else if (flag == 1) {
	s << this->getTag() << endln;
    } 
}

Response*
ZeroLengthImpact3D::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
    return new ElementResponse(this, 1, resid);
  
  // tangent stiffness matrix
  else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0)
    return new ElementResponse(this, 2, stiff);
  
  else 
    return Element::setResponse(argv,argc,output);
}

int 
ZeroLengthImpact3D::getResponse(int responseID, Information &eleInfo)
{
    if (responseID == 1)
	 return eleInfo.setVector(this->getResistingForce());
    else if (responseID == 2)
	 return eleInfo.setMatrix(this->getTangentStiff());
    else
      return Element::getResponse(responseID, eleInfo);
}

// Private methods
// determine the secondary/primary pair in contact, and setup Vectors (N,T1,T2)
int
ZeroLengthImpact3D::contactDetect(void)
{
  			  	

	  
	  int transientgap; 
	  transientgap = 1;   // 1: transient gap; 0: dynamic gap

	  Vector  secondaryNd;
	  Vector  primaryNd;

      //+--------------+-----------------+----------------+----------------+---------------+
      // NOTES: some methods to get displacements from nodes
      //+--------------+-----------------+----------------+----------------+---------------+
      // getDisp() :         get commit(k-1) disp, will be commit(k) after commit
      // getTrialDisp():     get Trial(k) disp
      // getIncrDisp():      get Trial(k)-Commit(k-1), will be 0 after commit
      // getIncrDeltaDisp(): get Trial(k)-Trial(k-1),  will be 0 after commit
      //+--------------+-----------------+----------------+----------------+---------------+

	  if (transientgap) 
	  {  ///////////// for transient gap //////////////////////////

		   secondaryNd = nodePointers[0]->getCrds() + nodePointers[0]->getTrialDisp();
           primaryNd= nodePointers[1]->getCrds() + nodePointers[1]->getTrialDisp();
	  }  else {
         ///////////// for dynamic gap ////////////////////////////
    	  secondaryNd = nodePointers[0]->getCrds() + nodePointers[0]->getIncrDisp();
          primaryNd= nodePointers[1]->getCrds() + nodePointers[1]->getIncrDisp();
	  }
      
      double Xs=secondaryNd(0)  - origin(0);
      double Ys=secondaryNd(1)  - origin(1);
	  double Zs=secondaryNd(2);
      double Rs=sqrt(Xs*Xs +Ys*Ys); 

      double Xp=primaryNd(0) - origin(0);
	  double Yp=primaryNd(1) - origin(1);
      double Zp=primaryNd(2);

	  double Rp=sqrt(Xp*Xp +Yp*Yp);

			

	  switch (directionID) {



         case 0:  // circular contact plane

	  

				if (transientgap) {

					gap = Rs-Rp - initGap;

				} else {

                   gap= gap_n + Rs - Rp - initGap; // dynamic gap

				}



				if (gap< 0.0) 

				{  // Not in contact

					return 0;

				} else 	{ // contact occur, setup contact vectors

			

					N(0)   =  -Xp/Rp ;

					N(1)   =  -Yp/Rp ;

					N(2)   =   0 ;

					N(3)   =   Xp/Rp ;

					N(4)   =   Yp/Rp ;

					N(5)   =   0 ;



					T1(0)  =   0;

					T1(1)  =   0;

					T1(2)  =   1;

					T1(3)  =   0;

					T1(4)  =   0;

					T1(5)  =  -1;



					T2(0)  =  -Yp/Rp ;

					T2(1)  =   Xp/Rp ;

					T2(2)  =   0 ;

					T2(3)  =   Yp/Rp ;

					T2(4)  =  -Xp/Rp ;

					T2(5)  =   0 ;



					return 1; 

				}

			



	 	case 1:   // normal of primary plane pointing to +X direction

				if (transientgap) {

					gap= Xp -Xs - initGap;             // transient gap

				} else {

                    gap= gap_n + Xp - Xs - initGap;    // dynamic gap

				}



				if (gap< 0.0)   

				{// Not in contact

					return 0;

				} else { 

				// contact occur, setup contact vectors

			

					N(0)   =   1;

					N(1)   =   0 ;

					N(2)   =   0 ;

					N(3)   =  -1 ;

					N(4)   =   0 ;

					N(5)   =   0 ;



					T1(0)  =   0;

					T1(1)  =   1;

					T1(2)  =   0;

					T1(3)  =   0;

					T1(4)  =  -1;

					T1(5)  =   0;



					T2(0)  =   0 ;

					T2(1)  =   0 ;

					T2(2)  =   1 ;

					T2(3)  =   0 ;

					T2(4)  =   0 ;

					T2(5)  =  -1 ;



					return 1; 

				}

			



		case 2:  // normal of primary plane pointing to +Y direction

				if (transientgap) {

					gap= Yp - Ys - initGap;            // transient gap

				} else {

					gap= gap_n + Yp - Ys - initGap;    // dynamic gap

				}



				if (gap< 0.0)  { // Not in contact

					return 0;

				} else	{ // contact occur, setup contact vectors

				

					N(0)   =   0;

					N(1)   =   1 ;

					N(2)   =   0 ;

					N(3)   =   0 ;

					N(4)   =  -1 ;

					N(5)   =   0 ;



					T1(0)  =   0;

					T1(1)  =   0;

					T1(2)  =   1;

					T1(3)  =   0;

					T1(4)  =   0;

					T1(5)  =  -1;



					T2(0)  =   1 ;

					T2(1)  =   0 ;

					T2(2)  =   0 ;

					T2(3)  =  -1 ;

					T2(4)  =   0 ;

					T2(5)  =   0 ;



					return 1; 

				}



		case 3:   // normal of primary plane pointing to +Z direction

			//          ___________ 

            //         |           |

			//         | secondary |  

			//         |___________| 

			//         |           |

			//         |  primary  |

            //         |           |

			//          -----------

			// 

				if (transientgap) {

					gap= Zp - Zs - initGap;         // transient gap

				} else {

					gap= gap_n + Zp - Zs - initGap; // dynamic gap

				}





				if (gap < 0.0)   // Not in contact

					return 0;

				else {	

					N(0)   =   0 ;

					N(1)   =   0 ;

					N(2)   =   1 ;

					N(3)   =   0 ;

					N(4)   =   0 ;

					N(5)   =  -1 ;



					T1(0)  =   1;

					T1(1)  =   0;

					T1(2)  =   0;

					T1(3)  =  -1;

					T1(4)  =   0;

					T1(5)  =   0;



					T2(0)  =   0 ;

					T2(1)  =   1 ;

					T2(2)  =   0 ;

					T2(3)  =   0 ;

					T2(4)  =  -1 ;

					T2(5)  =   0 ;



 				return 1;

				}



	  default:

            opserr << "ERROR!!!! ZeroLengthContact3D::ZeroLengthContact3D - the only available contact directions are 0,1,2,3\n";

            return -1;

	}  // end switch directionID

 }
 
void
ZeroLengthImpact3D::formResidAndTangent( int tang_flag ) 
{



	// trial displacement vectors

 	Vector DispTrialS(3); // trial disp for secondary node

	Vector DispTrialP(3); // trial disp for primary node

	// trial frictional force vectors (in local coordinate)

    Vector t_trial(2);

    double TtrNorm;



    // Coulomb friction law surface

	double Phi;     



    int i, j;



    //zero stiffness and residual 

    stiff.Zero( ) ;

    resid.Zero( ) ;

   

	// detect contact and set flag

    ContactFlag = contactDetect(); 



	//opserr<<this->getTag()<< " ZeroLengthContact3D::ContactFlag=" << ContactFlag<<endln;





	if (ContactFlag == 1) // contacted

	{  

       // contact presure;

		//pressure = Kn*gap;   // Kn : normal penalty

		KnANDpressure(); 

  

		DispTrialS=nodePointers[0]->getTrialDisp();

        DispTrialP=nodePointers[1]->getTrialDisp();



       //nodal displacements 

        double ul[6];



		ul[0]=DispTrialS(0);

		ul[1]=DispTrialS(1);

		ul[2]=DispTrialS(2);

		ul[3]=DispTrialP(0);

		ul[4]=DispTrialP(1);

		ul[5]=DispTrialP(2);



		t_trial.Zero();

        xi.Zero();	



		for (i=0; i<6; i++){

			xi(0) += T1(i)*ul[i];

			xi(1) += T2(i)*ul[i];

		}



		// Compute trial shear force

		for (i=0; i<2; i++)  t_trial(i)=Kt * (xi(i)-stickPt(i));  //Kt: tangential penalty

        TtrNorm=t_trial.Norm();



        // Coulomb friction law, trial state

		Phi = TtrNorm - (fs * pressure + cohesion);   // add cohesion

		

		if (Phi <= 0 ) { // stick case

  		//opserr<< "stick ...." << endln;



 			if ( tang_flag == 1 ) {

		    // stiff

				for (i=0; i<6; i++) {

					for (j=0; j<6; j++) {

						stiff(i,j) = Kn*(N(i)*N(j)) + Kt*(T1(i)*T1(j)+T2(i)*T2(j));	

					}

				}

			} //endif tang_flag

			// force

			for (i=0; i<6; i++) 

				 resid(i)= (-1*pressure)*N(i) + t_trial(0)*T1(i) + t_trial(1)*T2(i) ;

			//	 resid(i)= (-1*pressure)*N(i) - t_trial(0)*T1(i) - t_trial(1)*T2(i) ;



		} // end if stick

		else {              // slide case, non-symmetric stiff

            ContactFlag=2;  // set the contactFlag for sliding



		//	opserr<< "sliding ...." << endln;



            if ( tang_flag == 1 ) {

				// stiff

				double Pt1, Pt2;

				Pt1=t_trial(0)/TtrNorm;

				Pt2=t_trial(1)/TtrNorm;

				double C1=fs*Kn;

				double C2=Kt*(fs*pressure+cohesion)/TtrNorm;  // add cohesion, sept. 7, 2005



				for (i=0; i<6; i++) {

					for (j=0; j<6; j++) {

						stiff(i,j) = Kn*(N(i)*N(j)) - C1*(Pt1*T1(i)*N(j)+Pt2*T2(i)*N(j))

							    + C2*( (1-Pt1*Pt1)*T1(i)*T1(j) -    Pt1*Pt2 *T1(i)*T2(j)

					            - Pt1*Pt2 *T2(i)*T1(j) + (1-Pt1*Pt2)*T2(i)*T2(j)  );

					} //endfor i

				} //endfor j

			} // endif tang_flag

		

			// force

			double t1, t2;

			t1 = (fs*pressure+cohesion) * t_trial(0)/TtrNorm;  // add cohesion

            t2 = (fs*pressure+cohesion) * t_trial(1)/TtrNorm;  // add cohesion



			//opserr<<"gap=" << gap <<endln;

			//opserr<<"pressure= "<<pressure <<endln;



			for (i=0; i<6; i++) {

				resid(i) = (-1*pressure)*N(i)+t1*T1(i)+t2*T2(i) ;

			//	resid(i) = (-1*pressure)*N(i)-t1*T1(i)-t2*T2(i) ;

			}

		} //endif slide



	}  // endif ContactFlag



	// for NOT contact, do nothing, stiff and resid are zeroes



}












void
ZeroLengthImpact3D::KnANDpressure(void)
{
	
	gapT = gap;
	gapD = gapT - gapC;

	if (gapT<=0.0){ // not in contact
		pressT = 0.0;
		tangentT = 0.0;
	}
	if (gapT>0.0){ // contacted

		if (gapD>0.0) { // loading
			pressT = pressC+(Kn1*gapD);
			tangentT = Kn1;
			if ( (pressC+(Kn1*gapD)) > ((Kn1*Delta_y)+Kn2*(gapT - Delta_y)) ) {
				// Tstress = (K1*Delta_y)+K2*(Tstrain-gap-Delta_y); // original statement in impact material code
				pressT = (Kn1*Delta_y) + Kn2*(gapT - Delta_y);
				// Ttangent = K2; // original statement in impact material code
				tangentT = Kn2;
			}
		}
		if (gapD<0.0){ // unloading
			pressT = pressC+(Kn1*gapD);
			tangentT = Kn1;
			if ( (pressC+(Kn1*gapD)) < (Kn2*(gapT)) ) {
				// Tstress = K2*(Tstrain-gap);  // original statement in impact material code
				pressT = Kn2*(gapT);
				// Ttangent = K2;  // original statement in impact material code
				tangentT = Kn2;
			}
		}

		
	}

	pressure = pressT;
	Kn = tangentT;


	//opserr<<this->getTag()<< " ZeroLengthImpact3D::gapC=" << gapC <<endln;
	//opserr<<this->getTag()<< " ZeroLengthImpact3D::gapT=" << gapT <<endln;
	//opserr<<this->getTag()<< " ZeroLengthImpact3D::gapD=" << gapD <<endln;
	//opserr<<this->getTag()<< " ZeroLengthImpact3D::gapD*Kn1=" << gapD*Kn1 <<endln;
	//opserr<<this->getTag()<< " ZeroLengthImpact3D::pressT=***************" << pressT <<endln;
	//opserr<<this->getTag()<< " ZeroLengthImpact3D::Kn*gap=***************" << Kn*gap <<endln;

	
//                   
//                         
//                  |<  initGap  >|<  Delta_y >|   
//                  |                              Kn2
//           press  |                          +----------+
//      (pressure)  |                         /          /
//                  |                        /          /
//                  |                       /          /
//                  |                      /          /
//                  |                     /          /
//                  |                    /          /
//                  |                   /          /
//                  |                  /          /
//                  |                 /Kn1       /Kn1
//                  |                /          /
//                  |               /   Kn2    /
//                  |              /----------/
//         ---------+-------------+----------+-------------------
//                  |             |-> start of gap
//                  |   

}
