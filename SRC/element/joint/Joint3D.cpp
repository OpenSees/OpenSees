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

// $Revision: 1.7 $
// $Date: 2010-04-23 22:53:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/Joint3D.cpp,v $

// Written: Arash Altoontash, Gregory Deierlein
// Created: 03/02
// Revision: Arash

// Joint3D.cpp: implementation of the Joint3D class.
//
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <MP_Constraint.h>
#include <MP_Joint3D.h>
#include <ElementResponse.h>
#include <UniaxialMaterial.h>
#include <Joint3D.h>
#include <elementAPI.h>


Matrix Joint3D::K(45, 45);
Vector Joint3D::V(45);


void * OPS_ADD_RUNTIME_VPV(OPS_Joint3D)
{
    if (OPS_GetNDM() != 3 || OPS_GetNDF() != 6) {
  opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with Joint3D element\n";
  return 0;
    }

    // check the number of arguments is correct
    if (OPS_GetNumRemainingInputArgs() != 12 && OPS_GetNumRemainingInputArgs() != 16 ) {
  opserr << "WARNING incorrect number of arguments\n";
  //printCommand(argc, argv);
  opserr << "Want:\n";
  opserr << "element Joint3D Tag? NodI? NodJ? NodK? NodL? NodM? NodN? NodC? MatX? MatY? MatZ? LrgDsp?\n";
  opserr << "or:\n";
  opserr << "element Joint3D Tag? NodI? NodJ? NodK? NodL? NodM? NodN? NodC? MatX? MatY? MatZ? LrgDsp? -damage DmgX DmgY DmgZ\n";
  return 0;
    }

    // get the id and end nodes
    int idata[8];
    int num = 8;
    if (OPS_GetIntInput(&num, idata) < 0) {
  opserr << "WARNING invalid Joint3D int inputs" << endln;
  return 0;
    }
    int Joint3DId = idata[0];
    int iNode = idata[1];
    int jNode = idata[2];
    int kNode = idata[3];
    int lNode = idata[4];
    int mNode = idata[5];
    int nNode = idata[6];;

    // Get the center node
    int CenterNodeTag = idata[7];

    // check domain for existence of internal node tag
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return 0;
    Node *CenterNode = theDomain->getNode(CenterNodeTag);
    if (CenterNode != 0) {
  opserr << "WARNING node tag specified for the center node already exists.\n";
  opserr << "Use a new node tag.\n";
  opserr << "Joint3D element: " << Joint3DId << endln;
  return 0;
    }

    UniaxialMaterial *MatX = NULL;
    int MatXid;
    num = 1;
    if (OPS_GetIntInput(&num, &MatXid) < 0) {
  opserr << "WARNING invalid material ID for spring X\n";
  opserr << "Joint3D element: " << Joint3DId << endln;
  return 0;
    }

    MatX = OPS_getUniaxialMaterial(MatXid);
    if ( MatX == NULL )
    {
  opserr << "WARNING material not found\n";
  opserr << "Material: " << MatXid;
  opserr << "\nJoint3D element: " << Joint3DId << endln;
  return 0;
    }

    UniaxialMaterial *MatY = NULL;
    int MatYid;
    num = 1;
    if (OPS_GetIntInput(&num, &MatYid) < 0) {
  opserr << "WARNING invalid material ID for spring Y\n";
  opserr << "Joint3D element: " << Joint3DId << endln;
  return 0;
    }

    MatY = OPS_getUniaxialMaterial(MatYid);
    if ( MatY == NULL )
    {
  opserr << "WARNING material not found\n";
  opserr << "Material: " << MatYid;
  opserr << "\nJoint3D element: " << Joint3DId << endln;
  return 0;
    }

    UniaxialMaterial *MatZ = NULL;
    int MatZid;
    num = 1;
    if (OPS_GetIntInput(&num, &MatZid) < 0) {
  opserr << "WARNING invalid material ID for spring Z\n";
  opserr << "Joint3D element: " << Joint3DId << endln;
  return 0;
    }

    MatZ = OPS_getUniaxialMaterial(MatZid);
    if ( MatZ == NULL )
    {
  opserr << "WARNING material not found\n";
  opserr << "Material: " << MatZid;
  opserr << "\nJoint3D element: " << Joint3DId << endln;
  return 0;
    }

    int LargeDisp;
    num = 1;
    if (OPS_GetIntInput(&num, &LargeDisp) < 0) {
  // use 0 as default
  LargeDisp = 0;
    }


    Joint3D *theJoint3D;
    // Decide to use which constructor, based on the number of arguments
    if (OPS_GetNumRemainingInputArgs() == 12 ) {

  // Using Joint3D constructor without damage
      UniaxialMaterial* springModels[3] = { MatX, MatY, MatZ };
  theJoint3D = new Joint3D( Joint3DId,
          iNode,jNode,kNode,lNode,mNode,nNode,CenterNodeTag,
          springModels, theDomain, LargeDisp);

  // if get here we have successfully created the element and added it to the domain
  return theJoint3D;
    }

    else 			// if ( (argc-argStart) == 16  )
    {
  opserr<< "WARNING Using Joint3D constructor with damage not implemented in this version\n";
  return 0;
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


Joint3D::Joint3D()
  :Element(0, ELE_TAG_Joint3D),
  ExternalNodes(7), InternalConstraints(6),
  TheDomain(0), numDof(0), nodeDbTag(0), dofDbTag(0)
{
  for (int i = 0; i < 7; i++)
    theNodes[i] = NULL;

  for (int j = 0; j < 3; j++)
    theSprings[j] = NULL;
}



Joint3D::Joint3D(int tag, int nd1, int nd2, int nd3, int nd4, int nd5, int nd6, int IntNodeTag,
  UniaxialMaterial* springModels[],
  Domain* theDomain, int LrgDisp)
  :Element(tag, ELE_TAG_Joint3D), ExternalNodes(7), InternalConstraints(6),
  TheDomain(0), numDof(0), nodeDbTag(0), dofDbTag(0)
{
  int i;
  numDof = 45;

  K.Zero();
  V.Zero();

  TheDomain = theDomain;
  if (TheDomain == NULL) {
    opserr << "WARNING Joint3D(): Specified domain does not exist , Domain = 0\n";
    return;
  }

  // Save external node id's
  ExternalNodes(0) = nd1;
  ExternalNodes(1) = nd2;
  ExternalNodes(2) = nd3;
  ExternalNodes(3) = nd4;
  ExternalNodes(4) = nd5;
  ExternalNodes(5) = nd6;
  ExternalNodes(6) = IntNodeTag;


  // get  the external nodes
  for (i = 0; i < 6; i++)
  {
    theNodes[i] = NULL;
    theNodes[i] = TheDomain->getNode(ExternalNodes(i));
    if (theNodes[i] == NULL) {
      opserr << "WARNING Joint3D::setDomain(): Nd" << (i + 1) << ": ";
      opserr << ExternalNodes(i) << "does not exist in model for element \n" << *this;
      return;
    }
  }

  // check for a two dimensional domain, since this element supports only two dimensions 
  const Vector& end1Crd = theNodes[0]->getCrds();
  const Vector& end2Crd = theNodes[1]->getCrds();
  const Vector& end3Crd = theNodes[2]->getCrds();
  const Vector& end4Crd = theNodes[3]->getCrds();
  const Vector& end5Crd = theNodes[4]->getCrds();
  const Vector& end6Crd = theNodes[5]->getCrds();

  int dimNd1 = end1Crd.Size();
  int dimNd2 = end2Crd.Size();
  int dimNd3 = end3Crd.Size();
  int dimNd4 = end4Crd.Size();
  int dimNd5 = end5Crd.Size();
  int dimNd6 = end6Crd.Size();

  if (dimNd1 != 3 || dimNd2 != 3 || dimNd3 != 3 || dimNd4 != 3 || dimNd5 != 3 || dimNd6 != 3) {
    opserr << "WARNING Joint3D::setDomain(): has incorrect space dimension \n";
    opserr << "                                    space dimension not supported by Joint3D";
    return;
  }

  // now verify the number of dof at node ends
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();
  int dofNd3 = theNodes[2]->getNumberDOF();
  int dofNd4 = theNodes[3]->getNumberDOF();
  int dofNd5 = theNodes[4]->getNumberDOF();
  int dofNd6 = theNodes[5]->getNumberDOF();

  if (dofNd1 != 6 || dofNd2 != 6 || dofNd3 != 6 || dofNd4 != 6 || dofNd5 != 6 || dofNd6 != 6) {
    opserr << "WARNING Joint3D::Joint3D: has incorrect degrees of freedom \n";
    opserr << "                                    DOF not supported by Joint3D";
    return;
  }

  // check the joint size. The joint size must be non-zero
  Vector Center1(end1Crd);
  Vector Center2(end3Crd);
  Vector Center3(end5Crd);
  Center1 = Center1 - end2Crd;
  Center2 = Center2 - end4Crd;
  Center3 = Center3 - end6Crd;

  double L1 = Center1.Norm();
  double L2 = Center2.Norm();
  double L3 = Center3.Norm();

  if (Center1.Norm() < 1e-12 || Center2.Norm() < 1e-12 || Center3.Norm() < 1e-12) {
    opserr << "WARNING Joint3D::(): zero length\n";
    return;
  }

  // check if nodes are not located on each other and they can construct
  // a parallelogram
  Center1 = end1Crd + end2Crd;
  Center2 = end3Crd + end4Crd;
  Center3 = end5Crd + end6Crd;

  Center1 = 0.5 * Center1;
  Center2 = 0.5 * Center2;
  Center3 = 0.5 * Center3;

  Vector CenterTemp(Center2);
  CenterTemp = CenterTemp - Center1;
  if (CenterTemp.Norm() > 1e-6) {
    opserr << "WARNING Joint3D::(): can not construct a shear block over external nodes\n";
    opserr << "check the coordinates\n";
    return;
  }

  CenterTemp = Center3 - Center1;
  if (CenterTemp.Norm() > 1e-6) {
    opserr << "WARNING Joint3D::(): can not construct a shear block over external nodes\n";
    opserr << "check the coordinates\n";
    return;
  }

  // Generate internal node and add it up to domain
  theNodes[6] = new Node(IntNodeTag, 9, Center1(0), Center1(1), Center1(2));
  if (theNodes[6] == NULL) {
    opserr << "Joint3D::Joint3D - Unable to generate new nodes , out of memory\n";
  }
  else {
    if (TheDomain->addNode(theNodes[6]) == false)		// add intenal nodes to domain
      opserr << "Joint3D::Joint3D - unable to add internal nodeto domain\n";
  }

  // make copy of the uniaxial materials for the element

  if (springModels[0] == NULL) {
    opserr << "ERROR Joint3D::Joint3D(): The rotational spring in y'z' plane does not exist ";
    exit(-1);
  }
  else { theSprings[0] = springModels[0]->getCopy(); }

  if (springModels[1] == NULL) {
    opserr << "ERROR Joint3D::Joint3D(): The rotational spring in x'z' plane does not exist ";
    exit(-1);
  }
  else { theSprings[1] = springModels[1]->getCopy(); }

  if (springModels[2] == NULL) {
    opserr << "ERROR Joint3D::Joint3D(): The rotational spring in x'y' plane does not exist ";
    exit(-1);
  }
  else { theSprings[2] = springModels[2]->getCopy(); }

  for (i = 0; i < 3; i++)
  {
    if (theSprings[i] == NULL) {
      opserr << "ERROR Joint3D::Joint3D(): Can not make copy of uniaxial materials, out of memory ";
      exit(-1);
    }
  }


  // Generate and add constraints to domain

  // create MP_Joint constraint node 1
  InternalConstraints(0) = addMP_Joint(TheDomain, ExternalNodes(6), ExternalNodes(0), ExternalNodes(5), 8, ExternalNodes(3), 7, LrgDisp);
  if (InternalConstraints(0) < 0) {
    opserr << "WARNING Joint3D::Joint3D(): can not generate ForJoint MP at node 1\n";
    return;
  }

  // create MP_Joint constraint node 2
  InternalConstraints(1) = addMP_Joint(TheDomain, ExternalNodes(6), ExternalNodes(1), ExternalNodes(5), 8, ExternalNodes(3), 7, LrgDisp);
  if (InternalConstraints(1) < 0) {
    opserr << "WARNING Joint3D::Joint3D(): can not generate ForJoint MP at node 2\n";
    return;
  }

  // create MP_Joint constraint node 3
  InternalConstraints(2) = addMP_Joint(TheDomain, ExternalNodes(6), ExternalNodes(2), ExternalNodes(1), 6, ExternalNodes(5), 8, LrgDisp);
  if (InternalConstraints(2) < 0) {
    opserr << "WARNING Joint3D::Joint3D(): can not generate ForJoint MP at node 3\n";
    return;
  }

  // create MP_Joint constraint node 4
  InternalConstraints(3) = addMP_Joint(TheDomain, ExternalNodes(6), ExternalNodes(3), ExternalNodes(1), 6, ExternalNodes(5), 8, LrgDisp);
  if (InternalConstraints(3) < 0) {
    opserr << "WARNING Joint3D::Joint3D(): can not generate ForJoint MP at node 4\n";
    return;
  }

  // create MP_Joint constraint node 5
  InternalConstraints(4) = addMP_Joint(TheDomain, ExternalNodes(6), ExternalNodes(4), ExternalNodes(3), 7, ExternalNodes(1), 6, LrgDisp);
  if (InternalConstraints(4) < 0) {
    opserr << "WARNING Joint3D::Joint3D(): can not generate ForJoint MP at node 3\n";
    return;
  }

  // create MP_Joint constraint node 6
  InternalConstraints(5) = addMP_Joint(TheDomain, ExternalNodes(6), ExternalNodes(5), ExternalNodes(3), 7, ExternalNodes(1), 6, LrgDisp);
  if (InternalConstraints(5) < 0) {
    opserr << "WARNING Joint3D::Joint3D(): can not generate ForJoint MP at node 3\n";
    return;
  }
}


Joint3D::~Joint3D()
{

  if (TheDomain != NULL)
  {
    MP_Constraint* Temp_MP;
    for (int i = 0; i < 6; i++)
    {
      Temp_MP = TheDomain->getMP_Constraint(InternalConstraints(i));

      if (Temp_MP != NULL)
      {
        TheDomain->removeMP_Constraint(InternalConstraints(i));
        delete Temp_MP;
      }
    }
    if (theNodes[6] != NULL)
    {
      int intnodetag = theNodes[6]->getTag();
      Node* theNode = TheDomain->removeNode(intnodetag);
      if (theNode != 0) // have to check against eurned node in case node already gone!
        delete theNode;
    }
  }

  for (int i = 0; i < 3; i++)
    if (theSprings[i] != NULL) delete theSprings[i];
}



void Joint3D::setDomain(Domain* theDomain)
{
  //Ckeck domain not null - invoked when object removed from a domain
  if (theDomain == 0) {
    for (int i = 0; i < 7; i++) theNodes[i] = NULL;
  }
  else {

    TheDomain = theDomain;
    this->DomainComponent::setDomain(theDomain);

    for (int i = 0; i < 7; i++)
      if (theNodes[i] == 0)  theNodes[i] = TheDomain->getNode(ExternalNodes(i));
  }
}


int Joint3D::addMP_Joint(Domain* theDomain, int RetNodeID, int ConNodeID, int RotNodeID, int Rdof, int DspNodeID, int Ddof, int LrgDispFlag)
{
  MP_Constraint* Temp_MP;

  // create MP_ForJoint constraint
  Temp_MP = new MP_Joint3D(theDomain, RetNodeID, ConNodeID, RotNodeID, Rdof, DspNodeID, Ddof, LrgDispFlag);

  if (Temp_MP == NULL)
  {
    opserr << "Joint3D::addMP_Joint - WARNING ran out of memory for MP_Joint3D MP_Constraint ";
    return -1;
  }
  // Add the multi-point constraint to the domain
  if (theDomain->addMP_Constraint(Temp_MP) == false)
  {
    opserr << "Joint3D::addMP_Joint - WARNING could not add equalDOF MP_Constraint to domain ";
    delete Temp_MP;
    return -2;
  }
  return Temp_MP->getTag();;
}

//////////////////////////////////////////////////////////////////////
// Public methods called, taken care of for 2D element subclasses
//////////////////////////////////////////////////////////////////////

int Joint3D::update(void)
{
  const Vector& dispC = theNodes[6]->getTrialDisp();

  int result = 0;

  for (int i = 0; i < 3; i++)
  {
    if (theSprings[i] != NULL) result = theSprings[i]->setTrialStrain(dispC(i + 6));
    if (result != 0) break;
  }

  return result;
}

int Joint3D::commitState()
{
  int result = 0;

  for (int i = 0; i < 3; i++)
  {
    if (theSprings[i] != NULL) result = theSprings[i]->commitState();
    if (result != 0) break;
  }

  return result;
}

int Joint3D::revertToLastCommit()
{
  int result = 0;

  for (int i = 0; i < 3; i++)
  {
    if (theSprings[i] != NULL) result = theSprings[i]->revertToLastCommit();
    if (result != 0) break;
  }

  return result;
}

int Joint3D::revertToStart(void)
{
  int result = 0;

  for (int i = 0; i < 3; i++)
  {
    if (theSprings[i] != NULL) result = theSprings[i]->revertToStart();
    if (result != 0) break;
  }

  return result;
}


int Joint3D::getNumExternalNodes(void) const
{
  return 7;
}

const ID& Joint3D::getExternalNodes(void)
{
  return ExternalNodes;
}

Node** Joint3D::getNodePtrs(void)
{
  return theNodes;
}

int Joint3D::getNumDOF(void)
{
  return numDof;
}

const Matrix& Joint3D::getTangentStiff(void)
{
  double Ktangent[3];
  for (int i = 0; i < 3; i++)
  {
    Ktangent[i] = 0;
    if (theSprings[i] != NULL) Ktangent[i] = theSprings[i]->getTangent();
  }

  K.Zero();

  K(42, 42) = Ktangent[0];
  K(43, 43) = Ktangent[1];
  K(44, 44) = Ktangent[2];

  return K;
}


const Matrix& Joint3D::getInitialStiff(void)
{
  double Kintial[3];
  for (int i = 0; i < 3; i++)
  {
    Kintial[i] = 0;
    if (theSprings[i] != NULL) Kintial[i] = theSprings[i]->getTangent();
  }

  K.Zero();

  K(42, 42) = Kintial[0];
  K(43, 43) = -Kintial[0];
  K(44, 44) = Kintial[1];

  return K;
}


const Matrix& Joint3D::getDamp(void)
{
  K.Zero();
  return K;
}

const Matrix& Joint3D::getMass(void)
{
  K.Zero();
  return K;
}

void Joint3D::Print(OPS_Stream& s, int flag)
{
  s << "\nElement: " << getTag() << " type: Joint3D iNode: "
    << ExternalNodes(0) << " jNode: " << ExternalNodes(1) << "\n"
    << " kNode: " << ExternalNodes(2) << " lNode: " << ExternalNodes(3) << "\n"
    << " mNode: " << ExternalNodes(4) << " nNode: " << ExternalNodes(5) << "\n"
    << " Internal node: " << ExternalNodes(6) << "\n";
}

/////////////////////////////////////////////////////////////////////
// methods for applying and returning loads
//////////////////////////////////////////////////////////////////////

void Joint3D::zeroLoad(void)
{

}

int Joint3D::addLoad(ElementalLoad* theLoad, double loadFactor)
{
  return 0;
}

int Joint3D::addInertiaLoadToUnbalance(const Vector& accel)
{
  return 0;
}



const Vector& Joint3D::getResistingForce()
{
  double Force[3];
  for (int i = 0; i < 3; i++)
  {
    Force[i] = 0;
    if (theSprings[i] != NULL) Force[i] = theSprings[i]->getStress();
  }

  V.Zero();

  V(42) = Force[0];
  V(43) = Force[1];
  V(44) = Force[2];

  return V;
}

const Vector&
Joint3D::getResistingForceIncInertia()
{
  return this->getResistingForce();
}


int Joint3D::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
  // first determine the four corner points of the element based on
  // the display factor (a measure of the distorted image)
  // store this information in 2 3d vectors v1 and v2

  static Vector v1(3);
  static Vector v2(3);
  static Vector v3(3);
  static Vector v4(3);
  static Vector v5(3);
  static Vector v6(3);
  
  theNodes[0]->getDisplayCrds(v1, fact, displayMode);
  theNodes[1]->getDisplayCrds(v2, fact, displayMode);
  theNodes[2]->getDisplayCrds(v3, fact, displayMode);
  theNodes[3]->getDisplayCrds(v4, fact, displayMode);
  theNodes[2]->getDisplayCrds(v5, fact, displayMode);
  theNodes[3]->getDisplayCrds(v6, fact, displayMode);

  // draw the center lines
  int dummy;
  dummy = theViewer.drawLine(v1, v2, 1.0, 1.0);
  dummy = theViewer.drawLine(v3, v4, 1.0, 1.0);
  dummy = theViewer.drawLine(v5, v6, 1.0, 1.0);

  // calculate the eight corners of the block
  Vector va(3);
  Vector vb(3);
  Vector vc(3);

  va = v2 - v1;
  vb = v4 - v3;
  vc = v6 - v5;

  Vector vbegin(3);
  Vector vend(3);
  vbegin = v1 + 0.5 * vb - 0.5 * vc;
  vend = vbegin + va;
  dummy = theViewer.drawLine(vbegin, vend, 1.0, 1.0);

  vbegin = vend;
  vend = vbegin + vb;
  dummy = theViewer.drawLine(vbegin, vend, 1.0, 1.0);

  vbegin = vend;
  vend = vbegin - va;
  dummy = theViewer.drawLine(vbegin, vend, 1.0, 1.0);

  vbegin = vend;
  vend = vbegin - vb;
  dummy = theViewer.drawLine(vbegin, vend, 1.0, 1.0);

  vbegin = v1 - 0.5 * vb - 0.5 * vc;
  vend = vbegin + va;
  dummy = theViewer.drawLine(vbegin, vend, 1.0, 1.0);

  vbegin = vend;
  vend = vbegin + vb;
  dummy = theViewer.drawLine(vbegin, vend, 1.0, 1.0);

  vbegin = vend;
  vend = vbegin - va;
  dummy = theViewer.drawLine(vbegin, vend, 1.0, 1.0);

  vbegin = vend;
  vend = vbegin - vb;
  dummy = theViewer.drawLine(vbegin, vend, 1.0, 1.0);

  vbegin = v1 + 0.5 * vb - 0.5 * vc;
  vend = vbegin - vb;
  dummy = theViewer.drawLine(vbegin, vend, 1.0, 1.0);

  vbegin = v1 + 0.5 * vb + 0.5 * vc;
  vend = vbegin - vb;
  dummy = theViewer.drawLine(vbegin, vend, 1.0, 1.0);

  vbegin = v2 + 0.5 * vb - 0.5 * vc;
  vend = vbegin - vb;
  dummy = theViewer.drawLine(vbegin, vend, 1.0, 1.0);

  vbegin = v2 + 0.5 * vb + 0.5 * vc;
  vend = vbegin - vb;
  dummy = theViewer.drawLine(vbegin, vend, 1.0, 1.0);

  return 0;

}


//most-probably requires to be overridden
Response* Joint3D::setResponse(const char** argv, int argc, OPS_Stream& output)
{
  //
  // we compare argv[0] for known response types for the Truss
  //

  if (strcmp(argv[0], "node") == 0 || strcmp(argv[0], "internalNode") == 0)
    return new ElementResponse(this, 1, Vector(9));

  else if (strcmp(argv[0], "size") == 0 || strcmp(argv[0], "jointSize") == 0)
    return new ElementResponse(this, 2, Vector(3));

  else if (strcmp(argv[0], "moment") == 0 || strcmp(argv[0], "moments") == 0
    || strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0)
    return new ElementResponse(this, 3, Vector(3));

  else if (strcmp(argv[0], "defo") == 0 || strcmp(argv[0], "deformations") == 0 ||
    strcmp(argv[0], "deformation") == 0)
    return new ElementResponse(this, 4, Vector(3));

  else if (strcmp(argv[0], "defoANDforce") == 0 || strcmp(argv[0], "deformationANDforce") == 0 ||
    strcmp(argv[0], "deformationsANDforces") == 0)
    return new ElementResponse(this, 5, Vector(6));

  else if (strcmp(argv[0], "stiff") == 0 || strcmp(argv[0], "stiffness") == 0)
    return new ElementResponse(this, 6, Matrix(45, 45));

  else if (strcmp(argv[0], "plasticRotation") == 0 || strcmp(argv[0], "plasticDeformation") == 0)
    return new ElementResponse(this, 7, Vector(3));

  else
    return 0;

}

int Joint3D::getResponse(int responseID, Information& eleInformation)
{
  switch (responseID) {
  case -1:
    return -1;

  case 1:
    if (eleInformation.theVector != 0)
    {
      const Vector& disp = theNodes[6]->getTrialDisp();
      for (int i = 0; i < 9; i++)
        (*(eleInformation.theVector))(i) = disp(i);
    }
    return 0;

  case 2:
    if (eleInformation.theVector != 0)
    {
      const Vector& node1Crd = theNodes[0]->getCrds();
      const Vector& node2Crd = theNodes[1]->getCrds();
      const Vector& node3Crd = theNodes[2]->getCrds();
      const Vector& node4Crd = theNodes[3]->getCrds();
      const Vector& node5Crd = theNodes[4]->getCrds();
      const Vector& node6Crd = theNodes[5]->getCrds();

      const Vector& node1Disp = theNodes[0]->getDisp();
      const Vector& node2Disp = theNodes[1]->getDisp();
      const Vector& node3Disp = theNodes[2]->getDisp();
      const Vector& node4Disp = theNodes[3]->getDisp();
      const Vector& node5Disp = theNodes[4]->getDisp();
      const Vector& node6Disp = theNodes[5]->getDisp();

      Vector v1(3);
      Vector v2(3);
      Vector v3(3);
      Vector v4(3);
      Vector v5(3);
      Vector v6(3);

      // calculate the current coordinates of four external nodes
      for (int i = 0; i < 3; i++)
      {
        v1(i) = node1Crd(i) + node1Disp(i);
        v2(i) = node2Crd(i) + node2Disp(i);
        v3(i) = node3Crd(i) + node3Disp(i);
        v4(i) = node4Crd(i) + node4Disp(i);
        v5(i) = node5Crd(i) + node5Disp(i);
        v6(i) = node6Crd(i) + node6Disp(i);
      }
      v2 = v2 - v1;
      v4 = v4 - v3;
      v6 = v6 - v5;

      v1(0) = sqrt(v2(0) * v2(0) + v2(1) * v2(1) + v2(2) * v2(2));
      v1(1) = sqrt(v4(0) * v4(0) + v4(1) * v4(1) + v4(2) * v4(2));
      v1(2) = sqrt(v6(0) * v6(0) + v6(1) * v6(1) + v6(2) * v6(2));

      *(eleInformation.theVector) = v1;
    }
    return 0;

  case 3:
    if (eleInformation.theVector != 0)
    {
      for (int i = 0; i < 3; i++)
      {
        (*(eleInformation.theVector))(i) = 0.0;
        if (theSprings[i] != NULL)
          (*(eleInformation.theVector))(i) = theSprings[i]->getStress();
      }
    }
    return 0;

  case 4:
    if (eleInformation.theVector != 0)
    {
      for (int i = 0; i < 3; i++)
      {
        (*(eleInformation.theVector))(i) = 0.0;
        if (theSprings[i] != NULL)
          (*(eleInformation.theVector))(i) = theSprings[i]->getStrain();
      }
    }
    return 0;

  case 5:
    if (eleInformation.theVector != 0)
    {
      for (int i = 0; i < 3; i++)
      {
        (*(eleInformation.theVector))(i) = 0.0;
        (*(eleInformation.theVector))(i + 3) = 0.0;
        if (theSprings[i] != NULL)
        {
          (*(eleInformation.theVector))(i) = theSprings[i]->getStrain();
          (*(eleInformation.theVector))(i + 3) = theSprings[i]->getStress();
        }
      }
    }
    return 0;

  case 6:
    return eleInformation.setMatrix(this->getTangentStiff());

  case 7:
    if (eleInformation.theVector != 0)
    {
      for (int i = 0; i < 3; i++)
      {
        (*(eleInformation.theVector))(i) = 0.0;
        if (theSprings[i] != NULL && theSprings[i]->getInitialTangent() != 0.0)
        {
          (*(eleInformation.theVector))(i) =
            theSprings[i]->getStrain() - theSprings[i]->getStress() / theSprings[i]->getInitialTangent();
        }

      }
    }
    return 0;

  default:
    return -1;
  }
}


int Joint3D::sendSelf(int commitTag, Channel& theChannel)
{
  return 0;
}

int Joint3D::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  return 0;
}

