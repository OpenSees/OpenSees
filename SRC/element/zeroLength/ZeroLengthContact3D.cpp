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
                                                                        
// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLengthContact3D.cpp,v $
// $Revision: 1.3 $
// $Date: 2008-06-13 21:01:04 $

                                                                        
// Written: Gang Wang  (wang@ce.berkeley.edu)
//          Prof. Nicholas Sitar (nsitar@ce.berkeley.edu)
//
// Created: 27/08/2003
//
// Description: This file contains the implementation for the ZeroLengthContact3D class.


#include "ZeroLengthContact3D.h"
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

 
//static data
const int ZeroLengthContact3D::numberNodes = 2 ;

// static data for 3D
Matrix  ZeroLengthContact3D::stiff(3*numberNodes,3*numberNodes) ;
Vector  ZeroLengthContact3D::resid(3*numberNodes) ;
Matrix  ZeroLengthContact3D::zeroMatrix(3*numberNodes,3*numberNodes) ;

//*********************************************************************
//  Full Constructor:

ZeroLengthContact3D::ZeroLengthContact3D(int tag,
					 int Nd1, int Nd2, 
					 int direction, double Knormal, double Ktangent, 
					 double frictionRatio, double c, double origX, double origY )
  :Element(tag,ELE_TAG_ZeroLengthContact3D),     
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
  origin(0) = origX;
  origin(1) = origY;
  
  // set stick point cords in LOCAL basis
  stickPt(0)= 0;
  stickPt(1)= 0;
  
  // initialized contact flag be zero
  ContactFlag=0;
  
  gap_n = 0 ; 
}



//null constructor
ZeroLengthContact3D::ZeroLengthContact3D(void)
  :Element(0,ELE_TAG_ZeroLengthContact3D),     
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


//  Destructor:
//  delete must be invoked on any objects created by the object
//  and on the matertial object.
ZeroLengthContact3D::~ZeroLengthContact3D()
{
  if (load != 0)
    delete load;
  
  if (Ki != 0)
    delete Ki;
}


int
ZeroLengthContact3D::getNumExternalNodes(void) const
{
return 2;
}


const ID &
ZeroLengthContact3D::getExternalNodes(void) 
{
  return connectedExternalNodes;
}



Node **
ZeroLengthContact3D::getNodePtrs(void) 
{
  return nodePointers;
}

int
ZeroLengthContact3D::getNumDOF(void) 
{
  return numDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the ZeroLengthContact3D element
void
ZeroLengthContact3D::setDomain(Domain *theDomain)
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
 	"dofs at nodes in " << dofNd1 << " d problem\n"; 
     return;
    }
}   	 




int
ZeroLengthContact3D::commitState()
{

  // need to update stick point here 
  if (ContactFlag == 2 )   // update stick point only for slide case
    stickPt=xi; 

  // update gap for "dynamic gap" method
    gap_n = gap; 

	return 0;
}

int
ZeroLengthContact3D::revertToLastCommit()
{

	///////////////////////////////////////////
    // need to revert the stickPoint??
	xi=stickPt;

	return 0;
}


int
ZeroLengthContact3D::revertToStart()
{   

	// need to rezero stickPoint??
	stickPt.Zero();  
	return 0;
}


// calculate stress-strain relation -- M. Frank
/*
int
ZeroLengthContact3D::update(void)
{
 	return 0;
}
*/

const Matrix &
ZeroLengthContact3D::getTangentStiff(void)
{

  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;

}


const Matrix &
ZeroLengthContact3D::getInitialStiff(void)
{
  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}
    

const Matrix &
ZeroLengthContact3D::getDamp(void)
{
    // no mass 
 	zeroMatrix.Zero(); 
	return zeroMatrix;
}


const Matrix &
ZeroLengthContact3D::getMass(void)
{
    // no mass 
 	zeroMatrix.Zero(); 
	return zeroMatrix;
}


void 
ZeroLengthContact3D::zeroLoad(void)
{
  // does nothing now
}

int 
ZeroLengthContact3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // meaningless to addLoad to a contact !
  return 0;
}

int 
ZeroLengthContact3D::addInertiaLoadToUnbalance(const Vector &accel)
{
  // does nothing as element has no mass yet!
  return 0;
}


const Vector &
ZeroLengthContact3D::getResistingForce()
{

  int tang_flag = 0 ; //don't get the tangent
  formResidAndTangent( tang_flag ) ;

  return resid ;   
}


const Vector &
ZeroLengthContact3D::getResistingForceIncInertia()
{	
  // there is no Inertia 
 
  int tang_flag = 0 ; //don't get the tangent
  formResidAndTangent( tang_flag ) ;

  return  resid ;   
}


int
ZeroLengthContact3D::sendSelf(int commitTag, Channel &theChannel)
{
   // doing nothing here
	return 0;
}

int
ZeroLengthContact3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
// doing nothing here
	return 0;
}


int
ZeroLengthContact3D::displaySelf(Renderer &theViewer, int displayMode, float fact)
{ // nothing to display
    return 0;
}


void
ZeroLengthContact3D::Print(OPS_Stream &s, int flag)
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
ZeroLengthContact3D::setResponse(const char **argv, int argc, Information &eleInformation)
{
     if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
     return new ElementResponse(this, 1, resid);

     // tangent stiffness matrix
     else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0)
     return new ElementResponse(this, 2, stiff);

 	 else 
		return 0;
}


int 
ZeroLengthContact3D::getResponse(int responseID, Information &eleInfo)
{
 if (responseID == 1)
	 return eleInfo.setVector(this->getResistingForce());
 else if (responseID == 2)
	 return eleInfo.setMatrix(this->getTangentStiff());
 else
	 return -1;
}


// Private methods
// determine the slave/master pair in contact, and setup Vectors (N,T1,T2)
 int ZeroLengthContact3D::contactDetect(void)
 {
  			  	

	  
	  int transientgap; 
	  transientgap = 1;   // 1: transient gap; 0: dynamic gap

	  Vector  slaveNd;
	  Vector  masterNd;

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

		   slaveNd = nodePointers[0]->getCrds() + nodePointers[0]->getTrialDisp();
           masterNd= nodePointers[1]->getCrds() + nodePointers[1]->getTrialDisp();
	  }  else {
         ///////////// for dynamic gap ////////////////////////////
    	  slaveNd = nodePointers[0]->getCrds() + nodePointers[0]->getIncrDisp();
          masterNd= nodePointers[1]->getCrds() + nodePointers[1]->getIncrDisp();
	  }
      
      double Xs=slaveNd(0)  - origin(0);
      double Ys=slaveNd(1)  - origin(1);
	  double Zs=slaveNd(2);
      double Rs=sqrt(Xs*Xs +Ys*Ys); 

      double Xm=masterNd(0) - origin(0);
	  double Ym=masterNd(1) - origin(1);
      double Zm=masterNd(2);

	  double Rm=sqrt(Xm*Xm +Ym*Ym);

			

	  switch (directionID) {



         case 0:  // circular contact plane

	  

				if (transientgap) {

					gap = Rs-Rm;

				} else {

                   gap= gap_n + Rs - Rm; // dynamic gap

				}



				if (gap< 0) 

				{  // Not in contact

					return 0;

				} else 	{ // contact occur, setup contact vectors

			

					N(0)   =  -Xm/Rm ;

					N(1)   =  -Ym/Rm ;

					N(2)   =   0 ;

					N(3)   =   Xm/Rm ;

					N(4)   =   Ym/Rm ;

					N(5)   =   0 ;



					T1(0)  =   0;

					T1(1)  =   0;

					T1(2)  =   1;

					T1(3)  =   0;

					T1(4)  =   0;

					T1(5)  =  -1;



					T2(0)  =  -Ym/Rm ;

					T2(1)  =   Xm/Rm ;

					T2(2)  =   0 ;

					T2(3)  =   Ym/Rm ;

					T2(4)  =  -Xm/Rm ;

					T2(5)  =   0 ;



					return 1; 

				}

			



	 	case 1:   // normal of master plane pointing to +X direction

				if (transientgap) {

					gap= Xm -Xs;             // transient gap

				} else {

                    gap= gap_n + Xm - Xs;    // dynamic gap

				}



				if (gap< 0)   

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

			



		case 2:  // normal of master plane pointing to +Y direction

				if (transientgap) {

					gap= Ym - Ys;            // transient gap

				} else {

					gap= gap_n + Ym - Ys;    // dynamic gap

				}



				if (gap<=0)  { // Not in contact

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



		case 3:   // normal of master plane pointing to +Z direction

			//          ___________ 

            //         |           |

			//         |   slave   |  

			//         |___________| 

			//         |           |

			//         |   Master  |

            //         |           |

			//          -----------

			// 

				if (transientgap) {

					gap= Zm - Zs;         // transient gap

				} else {

					gap= gap_n + Zm - Zs; // dynamic gap

				}





				if (gap < 0)   // Not in contact

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





void  ZeroLengthContact3D::formResidAndTangent( int tang_flag ) 

{



	// trial displacement vectors

 	Vector DispTrialS(3); // trial disp for slave node

	Vector DispTrialM(3); // trial disp for master node

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

		pressure = Kn*gap;   // Kn : normal penalty

  

		DispTrialS=nodePointers[0]->getTrialDisp();

        DispTrialM=nodePointers[1]->getTrialDisp();



       //nodal displacements 

        double ul[6];



		ul[0]=DispTrialS(0);

		ul[1]=DispTrialS(1);

		ul[2]=DispTrialS(2);

		ul[3]=DispTrialM(0);

		ul[4]=DispTrialM(1);

		ul[5]=DispTrialM(2);



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





