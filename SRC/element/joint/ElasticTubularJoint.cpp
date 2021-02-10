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
                                                                        
// $Revision: 1.5 $
// $Date: 2010-05-05 23:09:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/ElasticTubularJoint.cpp,v $
                                                                        
// Written: Kia & Alanjari
//
// Description: This file contains the implementation for the ElasticTubularJoint class.
//
// What: "@(#) ElasticTubularJoint.C, revA"


// we specify what header files we need
#include "ElasticTubularJoint.h"
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

#define ELE_TAG_ElasticTubularJoint 2519

static const int NUM_NODE =2;
static const int NUM_DOF  =6;

static int numElasticTubularJoint = 0;


void *
OPS_ElasticTubularJoint(void)
{

  if (numElasticTubularJoint == 0) {
    numElasticTubularJoint++;
    opserr<<"ElasticTubularJoint element - Written by Kia & Alanjari\n";
  }

  // get the id and end nodes 
  int    iData[3];
  double  dData[6];
  int    numData;
  
  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[0]) != 0) {
    opserr << "\n WARNING invalid ElasticTubularJoint Tag" << endln;
    return 0;
  }
  
  
  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
    opserr << "\n WARNING invalid iNode for ElasticTubularJoint " << iData[0] << endln;
    return 0;
  }
 
 numData = 1;
 if (OPS_GetIntInput(&numData, &iData[2]) != 0) {
   opserr << "\n WARNING invalid jNode for ElasticTubularJoint " << iData[0] << endln;
   return 0;
 }
 
 numData = 1;
 if (OPS_GetDoubleInput(&numData, &dData[0]) != 0) {
   opserr << "\n WARNING invalid  brace diameter for ElasticTubularJoint " << iData[0] << endln;
   return 0;
 }
 
 numData = 1;
 if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
   opserr  << "\n WARNING invalid  brace_angle for ElasticTubularJoint " << iData[0] << endln ;
   return 0;
 }
 
 numData = 1;
 if (OPS_GetDoubleInput(&numData, &dData[2]) != 0) {
   opserr << "\n WARNING invalid E  for ElasticTubularJoint " << iData[0] << endln ;
   return 0;
 }
 
 numData = 1;
 if (OPS_GetDoubleInput(&numData, &dData[3]) != 0) {
   opserr << "\n WARNING invalid  chord diameter for ElasticTubularJoint "<< iData[0] << endln ;
   return 0;
 }
 
 numData = 1;
 if (OPS_GetDoubleInput(&numData, &dData[4]) != 0) {
   opserr  << "\n WARNING invalid  chord thickness for ElasticTubularJoint " <<iData[0] << endln ;
   return 0;
 }
 
 numData = 1;
 if (OPS_GetDoubleInput(&numData, &dData[5]) != 0) {
   opserr  << "\n WARNING invalid  chord angle for ElasticTubularJoint " <<iData[0] << endln ;
   return 0;
 }
 
 // now create the truss and add it to the Domain
 // Element *theTubularJoint = new ElasticTubularJoint(iData[0], iData[1], iData[2],dData[3] ,dData[2] , dData[5] , dData[0] ,dData[1] ,dData[4]);
 
 Element *theTubularJoint = new ElasticTubularJoint(iData[0], iData[1], iData[2],dData[0] ,dData[1] , dData[2] , dData[3] ,dData[4] ,dData[5]);
 
 if (theTubularJoint == 0) {
   opserr << "WARNING ran out of memory creating element with tag " << iData[0] << endln;
   return 0;
 }
 
 return theTubularJoint;
}


// typical constructor
ElasticTubularJoint::ElasticTubularJoint(int tag,int iNode , int jNode,
					 double Brace_Diameter,double Brace_Angle,
					 double e , double Chord_Diameter , double Chord_Thickness , double Chord_Angle )
 :Element(tag,ELE_TAG_ElasticTubularJoint), 
  l(0) ,cs(0) , sn(0), E(e), 
  braceD(Brace_Diameter), braceangle(Brace_Angle),
  chordD(Chord_Diameter),chordT(Chord_Thickness),chordangle(Chord_Angle),
  k(6,6),p(6),displacement(6),
  connectedExternalNodes(NUM_NODE)
{	
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2) 
    {
      opserr << "FATAl ElasticTubularJoint::ElasticTubularJoint - " <<  tag << "failed to create an ID of size 2\n";
      exit(-1);
    }
  connectedExternalNodes(0)=iNode;
  connectedExternalNodes(1)=jNode;
  
  // set node pointers to NUll
  theNodes[0]=0;
  theNodes[1]=0;
}

// constructor which should be invoked by an FE_ObjectBroker only
ElasticTubularJoint::ElasticTubularJoint()
:Element(0,ELE_TAG_ElasticTubularJoint), 
 l(0) ,cs(0) , sn(0),E(0), 
 braceD(0), braceangle(0),
 chordD(0),chordT(0),chordangle(0),
 k(6,6),p(6),displacement(6),
 connectedExternalNodes(NUM_NODE)
{   
  connectedExternalNodes(0)=0;
  connectedExternalNodes(1)=0;
  theNodes[0]=0;
  theNodes[1]= 0;
}

//  destructor - provided to clean up any memory
ElasticTubularJoint::~ElasticTubularJoint()
{
  
}

int
ElasticTubularJoint::getNumExternalNodes(void) const
{
  return NUM_NODE;
}

const ID &
ElasticTubularJoint::getExternalNodes(void) 
{
  return this->connectedExternalNodes;
}

Node **
ElasticTubularJoint::getNodePtrs(void) 
{
  return theNodes;
}

int
ElasticTubularJoint::getNumDOF(void) {
  return NUM_DOF ;
}



// method: setDomain()
//    to set a link to the enclosing Domain, ensure nodes exist in Domain
//    and set pointers to these nodes, also determines the length 

void
ElasticTubularJoint::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    return;
  }
  
  // first set the node pointers
  theNodes[0]= theDomain->getNode(connectedExternalNodes(0));
  theNodes[1]= theDomain->getNode(connectedExternalNodes(1));
  if (theNodes[0]==0 )
    {
      opserr << "  Node " << this->connectedExternalNodes(0) << " does not exit in the domain" << endln;
      return ; // don't go any further
    }
  
  if (theNodes[1]==0) 
    {
      opserr << "  Node " << this->connectedExternalNodes(1) << " does not exit in the domain  " << endln;
      
      return ; // don't go any further
    }
  
  // now determine the number of dof
  int nodf1= theNodes[0]->getNumberDOF();
  int nodf2= theNodes[1]->getNumberDOF();
  if (nodf1 != 3 || nodf2 != 3 )
    {
      opserr << "  3 dof required at each nodes " << endln;
      
      return ; // don't go any further
    }
  // call the base class method
  
  this->DomainComponent::setDomain(theDomain);
  
  // now determine the length
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();
  double dx = end2Crd(0)-end1Crd(0);
  double dy = end2Crd(1)-end1Crd(1);
  
  l = sqrt(dx*dx + dy*dy);
  if (l==0)
    {
      opserr << " \n ElasticTubularJoint " << this->getTag() << " has zero length" << endln;
      
      return ; // don't go any further
      
    }
  
  cs = dx/l;
  sn= dy/l;
  braceangle = braceangle *3.1415926535897932384626433832795/180;
  chordangle= (90 - chordangle)* 3.1415926535897932384626433832795/180;
}





int
ElasticTubularJoint::commitState()
{
  return 0;
}

int
ElasticTubularJoint::revertToLastCommit()
{
  return 0;
}

int
ElasticTubularJoint::revertToStart()
{
  return 0;
}


int
ElasticTubularJoint::update()
  
{
  return 0;
}




const Matrix &
ElasticTubularJoint::getTangentStiff(void)
{	
	
  if (l == 0.0)
    { // length = zero - problem in setDomain() warning message already K.Zero();	
      this->k.Zero();
      return  k;
    }
  
  double gama = this->chordD/(2* this->chordT);
  double betta= this->braceD/this->chordD;
  
  double res1= pow(gama ,2.15);
  double res2= pow((1-betta),1.3);
  double res3=pow((sin(braceangle)),0.19);
  TangLJFv =(1.95*res1*res2*res3/(this->E*this->chordD));
  
  double res4= pow(gama,1.73);
  double res5= exp(-4.52*betta);
  double res6= pow(sin(braceangle),1.22);
  TangLJFipb=(134*res4*res5*res6/(this->E*chordD*chordD*chordD));
  
  double a = TangLJFv;
  double c1 = TangLJFipb; 
  double b= 1/tan(3.1415926535897932384626433832795/2);
  
  double s = sin(chordangle);
  double c= cos(chordangle);
  double s2= pow(s , 2);
  double c2= pow(c ,2);
  
  k(0,0) = s2/b + c2/a;
  k(0,1)=k(1,0) = (a-b)/(a*b) *s*c;
  k(0,2)=k(2,0) = -(c2/a+s2/b)*l*sn +(a-b)/(a*b)*s*c*l*cs;
  k(0,3)= k(3,0) =-(s2/b+c2/a);
  k(0,4)= k(4,0)= -(a-b)/(a*b)*s*c;
  k(0,5) = k( 5,0) = 0;
  
  k(1,1) = c2/b + s2/a;
  k(1,2)= k(2,1)=(c2/b + s2/a)*l*cs - ((a-b)/(a*b)*s*c*l*sn);
  k(1,3)= k(3,1) = -(a-b)/(a*b)*s*c ;
  k(1,4) = k(4,1) = -(c2/b + s2/a );
  k(1,5) =k(5,1) =0 ;
  
  k(2,2)=((s2/b+c2/a)*l*sn-(a-b)/(a*b)*s*c*l*cs)*l*sn + ((c2/b+s2/a)*l*cs - (a-b)/(a*b)*s*c*l*sn)*l*cs + 1/c1;
  k(2,3) = k(3,2) = (s2/b + c2/a)*l*sn - (a-b)/(a*b)*s*c*l*cs ;
  k(2,4)= k(4,2) = (a-b)/(b*a)*s*c*l*sn - (c2/b + s2/a)*l*cs ;
  k(2,5) =k(5,2) = -1/c1;
  
  k(3,3) = s2/b + c2/a;
  k(3,4)=k(4,3) = (a-b)/(a*b)*s*c;
  k(3,5)=k(5,3)=0;
  
  k(4,4)= c2/b + s2/a ;
  k(4,5)=k(5,4)=0;
  k(5,5)= 1/c1;
  
  return k;
  
  
}


const Matrix &
ElasticTubularJoint::getInitialStiff(void)
{

  
  if ( l== 0) 
    {
      this->k.Zero();
      return k;
    }
  
  double gama = this->chordD/(2* this->chordT);
  double betta= this->braceD/this->chordD;
  
  double res1= pow(gama ,2.15);
  double res2= pow((1-betta),1.3);
  double res3=pow((sin(braceangle)),0.19);
  TangLJFv =(1.95*res1*res2*res3/(this->E*this->chordD));
  
  double res4= pow(gama,1.73);
  double res5= exp(-4.52*betta);
  double res6= pow(sin(braceangle),1.22);
  TangLJFipb=(134*res4*res5*res6/(this->E*chordD*chordD*chordD));
  
  double a = TangLJFv;
  double c1 = TangLJFipb; 
  double b= 1/tan(3.1415926535897932384626433832795/2);
  
  double s = sin(chordangle);
  double c= cos(chordangle);
  double s2= pow(s , 2);
  double c2= pow(c ,2);
  
  k(0,0) = s2/b + c2/a;
  k(0,1)=k(1,0) = (a-b)/(a*b) *s*c;
  k(0,2)=k(2,0) = -(c2/a+s2/b)*l*sn +(a-b)/(a*b)*s*c*l*cs;
  k(0,3)= k(3,0) =-(s2/b+c2/a);
  k(0,4)= k(4,0)= -(a-b)/(a*b)*s*c;
  k(0,5) = k( 5,0) = 0;
  
  k(1,1) = c2/b + s2/a;
  k(1,2)= k(2,1)=(c2/b + s2/a)*l*cs - ((a-b)/(a*b)*s*c*l*sn);
  k(1,3)= k(3,1) = -(a-b)/(a*b)*s*c ;
  k(1,4) = k(4,1) = -(c2/b + s2/a );
  k(1,5) =k(5,1) =0 ;
  
  k(2,2)=((s2/b+c2/a)*l*sn-(a-b)/(a*b)*s*c*l*cs)*l*sn + ((c2/b+s2/a)*l*cs - (a-b)/(a*b)*s*c*l*sn)*l*cs + 1/c1;
  k(2,3) = k(3,2) = (s2/b + c2/a)*l*sn - (a-b)/(a*b)*s*c*l*cs ;
  k(2,4)= k(4,2) = (a-b)/(b*a)*s*c*l*sn - (c2/b + s2/a)*l*cs ;
  k(2,5) =k(5,2) = -1/c1;
  
  k(3,3) = s2/b + c2/a;
  k(3,4)=k(4,3) = (a-b)/(a*b)*s*c;
  k(3,5)=k(5,3)=0;
  
  k(4,4)= c2/b + s2/a ;
  k(4,5)=k(5,4)=0;
  k(5,5)= 1/c1;
  
  return k;
}



const Vector &
ElasticTubularJoint::getResistingForce()
{
	
    
  const Vector &disp1 = theNodes[0]->getTrialDisp();
  const Vector &disp2 = theNodes[1]->getTrialDisp();
  
  displacement(0)=disp1(0);
  displacement(1)=disp1(1);
  displacement(2)=disp1(2);
  displacement(3)=disp2(0);
  displacement(4)=disp2(1);
  displacement(5)=disp2(2);
  
  p= k*displacement; 
  
  
  return p;
}
int 
ElasticTubularJoint::sendSelf(int Tag, Channel &theChannel)
{
  int res=0;
  static  Vector data(9);
  data(0)=this->getTag();
  data(1)= this->connectedExternalNodes(0);
  data(2)= this->connectedExternalNodes(1);
  data(3)= this->braceD;
  data(4)= this->braceangle;
  data(5)=this->E;
  data(6)=this->chordD;
  data(7)= this->chordT;
  data(8) = this->chordangle;
  // send the data vector
  res+=theChannel.sendVector(this->getDbTag(), Tag , data);
  if (res < 0 )
    {
      opserr << " ElasticTubularJoint::sendSlef--could not send data vector \n ";
      return res;
    }
  
  
  return res;
}

int 
ElasticTubularJoint::recvSelf(int Tag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res=0;
  static Vector data(9);
  res += theChannel.recvVector(this->getDbTag(), Tag, data);
  if (res < 0) 
    {
      opserr << " Tubular Joint Element ::recvself-- could not receive data vector \n ";
      return res ;
    }
  this->setTag((int)data(0));
  this->connectedExternalNodes(0)=(int)data(1);
  this->connectedExternalNodes(1)= (int)data(2);
  this->braceD=data(3);
  this->braceangle=data(4);
  this->E=data(5);
  this->chordD = data(6);
  this->chordT = data(7);
  this ->chordangle=data(8);
  
  return res;
}






void
ElasticTubularJoint::Print(OPS_Stream &s, int flag)
{
  s << " Element tag:" << this->getTag() << endln;
  s << "  iNode : " << this->connectedExternalNodes(0) << endln;
  s << "  jNode : " << this->connectedExternalNodes(1) << endln ;
  s << "  E : " << this->E << endln;
  s << "   Axial Stiffness =" <<1/( this->TangLJFv*pow(sin(braceangle),2));
  s << " In Plane Bending Stiffness = "  << 1/(this->TangLJFipb) << endln;
  
  s << " End 1 Forces (P,V,M) : " << " (" << p(0) << " , " << p(1) << " , " << p(2) << " )" << endln;
  
  s << " End 2 Forces (P,V,M) :" << " (" << p(3) << " ," << p(4) << " ," << p(5) << " )" << endln;
}






int ElasticTubularJoint::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  // ensure setDomain() worked
  if (l == 0.0)
    return 0;
  
  static Vector v1(3);
  static Vector v2(3);

  theNodes[0]->getDisplayCrds(v1, fact, displayMode);
  theNodes[1]->getDisplayCrds(v2, fact, displayMode);

  return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

