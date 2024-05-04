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

// Written by: Long Chen, Pedro Arduino (parduino@uw.edu), Wenyang Zhang and fmk
//
// Four node PML2D_12 element .. a c++ wrapper to fortran routine 
// provided by Wenyang Zhang (zwyll@ucla.edu), University of California, Los Angeles
//
// University of Washington, UC. Los Angeles, U.C. Berkeley, 12, 2020


#include "PML2D_12.h"

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <OPS_Globals.h>
#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Parameter.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>


void* OPS_PML2D_12()
{
  if (OPS_GetNumRemainingInputArgs() < (6+PML2D_12_NUM_PROPS)) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element PML2D_12 eleTag? [PML2D_12_NUM_NODES integer nodeTags] [PML2D_12_NUM_PROPS material properties]\n";
    return 0;
  }
  
  int idata[6];
  int num = 6;
  if (OPS_GetIntInput(&num,idata)<0) {
    opserr<<"WARNING: invalid integer data\n";
    return 0;
  }
  
  double dData[PML2D_12_NUM_PROPS]; num = PML2D_12_NUM_PROPS;
  // make last two indicies zero 
  dData[PML2D_12_NUM_PROPS-2] = 0.0;
  dData[PML2D_12_NUM_PROPS-1] = 0.0;

  if (OPS_GetDoubleInput(&num,dData) < 0) {
    opserr<<"WARNING: invalid double data\n";
    return 0;
  }	

  return new PML2D_12(idata[0],&idata[1],dData);
}

//static data
Matrix PML2D_12::tangent;
Vector PML2D_12::resid(PML2D_12_NUM_DOF);
int    PML2D_12::eleCount=0;    


//null constructor
PML2D_12::PML2D_12( ) 
:Element( 0, ELE_TAG_PML2D_12),
 connectedExternalNodes(PML2D_12_NUM_NODES)
{
  for (int i=0; i<PML2D_12_NUM_NODES; i++ ) {
    nodePointers[i] = 0;
  }
}


//*********************************************************************
//full constructor
PML2D_12::PML2D_12(int tag, 
	     int *nodeTags,
	     double *eleData)
  :Element(tag, ELE_TAG_PML2D_12),
   connectedExternalNodes(PML2D_12_NUM_NODES)
{
  eleCount++;
  if (eleCount == 1) {
    opserr << "Perfectly Matched Layer 2D_12 (PML) element -  Written: W. Zhang, E. Taciroglu, A. Pakzad, P. Arduino, UCLA, UCLA, U.Washington, U.Washington\n ";
  }
  for (int i=0; i<PML2D_12_NUM_NODES; i++) {
    connectedExternalNodes(i) = nodeTags[i];
    nodePointers[i] = 0;
  } 


  for (int i=0; i<PML2D_12_NUM_PROPS; i++)
    props[i] = eleData[i];

}

//destructor 
PML2D_12::~PML2D_12( )
{

}


//set domain
void  PML2D_12::setDomain( Domain *theDomain ) 
{  

  int i ;
 
  //node pointers
  for ( i=0; i<PML2D_12_NUM_NODES; i++ ) 
     nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;

  this->DomainComponent::setDomain(theDomain);


  // 
  // set constant matrices by invoking fortran routine
  //

  //  double coords[PML2D_12_NUM_NODES][2];
  double coords[PML2D_12_NUM_NODES*2];

  for (int i = 0; i < PML2D_12_NUM_NODES; i++ )  {
      const Vector &loc = nodePointers[i]->getCrds();
      coords[i*2] = loc(0);      
      coords[i*2+1] = loc(1);      
  } 

  int NDOFEL = PML2D_12_NUM_DOF;
  int NPROPS = PML2D_12_NUM_PROPS;
  int MCRD = 2; 
  int NNODE = 4;
  double G[PML2D_12_NUM_DOF * PML2D_12_NUM_DOF];

  double tempK[PML2D_12_NUM_DOF*PML2D_12_NUM_DOF];
  double tempC[PML2D_12_NUM_DOF*PML2D_12_NUM_DOF];
  double tempM[PML2D_12_NUM_DOF*PML2D_12_NUM_DOF];
  pml2d_(tempK, 
       tempC, 
       tempM,  
       G, 
       &NDOFEL, 
       props, 
       &NPROPS, 
       coords, 
       &MCRD,
       &NNODE);

  for (int i = 0; i < 20; i++) {
    for (int j = 0; j < 20; j++) {
      int index1 = i*20 + j;
      int ii, jj;
      if (i%5 <2) {
        ii = int(i/5)*2 + i%5;
      } else {
        ii = 8 + int(i/5)*3 + i%5 -2;
      }
      if (j%5 <2) {
        jj = int(j/5)*2 + j%5;
      } else {
        jj = 8 + int(j/5)*3 + j%5 -2;
      }

      int index2 = ii*20 + jj;
      K[index2] = tempK[index1];
      C[index2] = tempC[index1];
      M[index2] = tempM[index1];
    }
  }


}


//get the number of external nodes
int  PML2D_12::getNumExternalNodes( ) const
{
  return PML2D_12_NUM_NODES ;
} 
 

//return connected external nodes
const ID&  PML2D_12::getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 

Node **  
PML2D_12::getNodePtrs(void) 
{
  return nodePointers ;
} 

//return number of dofs
int  PML2D_12::getNumDOF( ) 
{
  return PML2D_12_NUM_DOF ;
}


//commit state
int  PML2D_12::commitState( )
{
  // call element commitState to do any base class stuff
  int success = 0;
  if ((success = this->Element::commitState()) != 0) {
    opserr << "PML2D_12::commitState () - failed in base class";
  }    

  return success ;
}
 


//revert to last commit 
int  PML2D_12::revertToLastCommit( ) 
{
  int success = 0;
  return success ;
}
    

//revert to start 
int  PML2D_12::revertToStart( ) 
{
  int success = 0 ;
  return success ;
}

//print out element data
void  PML2D_12::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {
    
    s << "#PML2D_12\n";
    
    int i;
    const int numNodes = PML2D_12_NUM_NODES;
    
    for (i = 0; i < numNodes; i++) {
      const Vector &nodeCrd = nodePointers[i]->getCrds();
      const Vector &nodeDisp = nodePointers[i]->getDisp();
      s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << nodeCrd(2)
	<< " " << nodeDisp(0) << " " << nodeDisp(1) << " " << nodeDisp(2) << endln;
    }
  }
  if (flag == OPS_PRINT_CURRENTSTATE) {
    
    s << "PML2D_12 \n";
    s << "Element Number: " << this->getTag() << endln;
    s << "Nodes: " << connectedExternalNodes;
    
    s << endln;
    s << this->getTag() << " " << connectedExternalNodes(0)
      << " " << connectedExternalNodes(1)
      << " " << connectedExternalNodes(2)
      << " " << connectedExternalNodes(3)
      << endln;
    
    s << "Resisting Force (no inertia): " << this->getResistingForce();
  }
    
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"PML2D_12\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
    for (int i = 1; i < 3; i++)
      s << connectedExternalNodes(i) << ", ";
    s << connectedExternalNodes(3) << "], ";
  }
}
 
 
//return stiffness matrix 
const Matrix&  PML2D_12::getTangentStiff( ) 
{
  tangent.setData(K, PML2D_12_NUM_DOF, PML2D_12_NUM_DOF);
  return tangent;
}    


const Matrix&  PML2D_12::getInitialStiff( ) 
{
  return this->getTangentStiff();
}    


//return mass matrix
const Matrix&  
PML2D_12::getMass( ) 
{
  tangent.setData(M, PML2D_12_NUM_DOF, PML2D_12_NUM_DOF);
  return tangent;
} 

//return damping matrix
const Matrix&  
PML2D_12::getDamp( ) 
{
  tangent.setData(C, PML2D_12_NUM_DOF, PML2D_12_NUM_DOF);
  return tangent;
} 



void  PML2D_12::zeroLoad( )
{
  return ;
}


int 
PML2D_12::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  return -1;
}

//get residual
const Vector&  
PML2D_12::getResistingForce( ) 
{
  int numNodalDOF = 2;
  static Vector theVector(PML2D_12_NUM_DOF);

  // get K into stiff
  tangent.setData(K, PML2D_12_NUM_DOF, PML2D_12_NUM_DOF);

  int loc = 0;
  for (int i = 0; i < 4; i++ )  {
    const Vector &uNode = nodePointers[i]->getTrialDisp();
    for (int j=0; j<numNodalDOF; j++)
      theVector(loc++) = uNode(j);      
  }
  const Vector &uNode = nodePointers[4]->getTrialDisp();
  for (int j=0; j<12; j++)
    theVector(loc++) = uNode(j);

  resid.addMatrixVector(0.0, tangent, theVector, 1.0);

  return resid ;   
}


//get residual
const Vector&  
PML2D_12::getResistingForceIncInertia( ) 
{

  int numNodalDOF = 2;
  static Vector theVector(PML2D_12_NUM_DOF);
  static Matrix theMatrix(PML2D_12_NUM_DOF, PML2D_12_NUM_DOF);

  //
  // perform: R = P(U) - Pext(t)
  //

  this->getResistingForce();

  //
  // perform: R = R - M * a
  //

  int loc = 0;
  Node **theNodes = this->getNodePtrs();
  for (int i=0; i<4; i++) {
    const Vector &acc = theNodes[i]->getTrialAccel();
    for (int j=0; j<numNodalDOF; j++) {
      theVector(loc++) = acc(j);
    }
  } 
  {
  const Vector &uNode = nodePointers[4]->getTrialAccel();
  for (int j=0; j<12; j++)
    theVector(loc++) = uNode(j);
  }
  resid.addMatrixVector(1.0, this->getMass(), theVector, +1.0);

  //
  // perform: R = R + (alphaM * M + betaK0 * K0 + betaK * K) * v
  //            = R + D * v
  //

  // determine the vel vector from ele nodes
  loc = 0;
  for (int i=0; i<4; i++) {
    const Vector &vel = theNodes[i]->getTrialVel();
    for (int j=0; j<numNodalDOF; j++) {
      theVector(loc++) = vel[j];
    }
  }
  {
  const Vector &uNode = nodePointers[4]->getTrialVel();
  for (int j=0; j<12; j++)
    theVector(loc++) = uNode(j);
  }

  // finally the C * v
  resid.addMatrixVector(1.0, this->getDamp(), theVector, 1.0);

  return resid;
}


int  
PML2D_12::update(void) 
{
  return 0;
}



int  PML2D_12::sendSelf (int commitTag, 
		      Channel &theChannel)
{
  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // PML2D_12 packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(PML2D_12_NUM_PROPS + 1);
  data(0) = this->getTag();

  for (int ii = 1; ii <= PML2D_12_NUM_PROPS; ii++) {
      data(ii) = props[ii - 1];
  }

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
      opserr << "WARNING PML2D_12::sendSelf() - " << this->getTag() << " failed to send Vector\n";
      return res;
  }


  // PML2D_12 then sends the tags of its four nodes
  res += theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
      opserr << "WARNING PML2D_12::sendSelf() - " << this->getTag() << " failed to send ID\n";
      return res;
  }

  return res;
  }
    
int  PML2D_12::recvSelf (int commitTag, 
		      Channel &theChannel, 
		      FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  // PML2D_12 creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(PML2D_12_NUM_PROPS + 1);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
      opserr << "WARNING PML2D_12::recvSelf() - failed to receive Vector\n";
      return res;
  }

  this->setTag((int)data(0));

  for (int ii = 1; ii <= PML2D_12_NUM_PROPS; ii++) {
      props[ii - 1] = data(ii);
  }

  // PML2D_12 now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
      opserr << "WARNING PML2D_12::recvSelf() - " << this->getTag() << " failed to receive ID\n";
      return res;
  }

  return res;
}


int
PML2D_12::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // get the end point display coords
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    static Vector v4(3);
    nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
    nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
    nodePointers[2]->getDisplayCrds(v3, fact, displayMode);
    nodePointers[3]->getDisplayCrds(v4, fact, displayMode);

    // place values in coords matrix
    static Matrix coords(4, 3);
    for (int i = 0; i < 3; i++) {
        coords(0, i) = v1(i);
        coords(1, i) = v2(i);
        coords(2, i) = v3(i);
        coords(3, i) = v4(i);
    }

    // fill RGB vector
    static Vector values(4);
    for (int i = 0; i < 4; i++)
        values(i) = 1.0;

    // draw the polygon
    return theViewer.drawPolygon(coords, values, this->getTag());
}

Response*
PML2D_12::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType","PML2D_12");
  output.attr("eleTag",this->getTag());
  for (int i=1; i<=8; i++) {
    sprintf(outputData,"node%d",i);
    output.attr(outputData,nodePointers[i-1]->getTag());
  }

  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {

    for (int i=1; i<=8; i++) {
      sprintf(outputData,"P1_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P2_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P3_%d",i);
      output.tag("ResponseType",outputData);
    }

    theResponse = new ElementResponse(this, 1, resid);
  }  
  output.endTag(); // ElementOutput
  return theResponse;
}

int 
PML2D_12::getResponse(int responseID, Information &eleInfo)
{
  static Vector stresses(48);

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  return -1;
}

int
PML2D_12::setParameter(const char **argv, int argc, Parameter &param)
{
  int res = -1;
  return res;
}
    
int
PML2D_12::updateParameter(int parameterID, Information &info)
{
    int res = -1;
    return res;
}

