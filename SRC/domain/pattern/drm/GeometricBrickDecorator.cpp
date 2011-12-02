#include <cstdlib>

#include "GeometricBrickDecorator.h"

GeometricBrickDecorator::GeometricBrickDecorator() 
{
  this->myBrick = NULL;
  this->myDomain = NULL;
}

GeometricBrickDecorator::~GeometricBrickDecorator()
{ }

void
GeometricBrickDecorator::setBrick(Element*  elem)
{
  this->myBrick = elem;
}

void
GeometricBrickDecorator::setDomain(Domain * theDomain)
{
  this->myDomain = theDomain;
}

void 
GeometricBrickDecorator::clearBrick()
{
  this->myBrick = NULL;
}

void 
GeometricBrickDecorator::clearDomain()
{
  this->myDomain = NULL;
}

void 
GeometricBrickDecorator::clearAll()
{
  this->clearBrick();
  this->clearDomain();
}

bool
GeometricBrickDecorator::isClear()
{
  return ( (this->myBrick == NULL) && (this->myDomain == NULL) );
}

bool
GeometricBrickDecorator::isEmpty()
{
  return ( (this->myBrick == NULL) || (this->myDomain == NULL) );
}

bool
GeometricBrickDecorator::isPointInVolume(const Vector &pP)
{
  if ( (this->isClear()) || (this->isEmpty()) ) {
    return false;
  }
  
  const ID extNodes = this->myBrick->getExternalNodes();

  int numExtNodes = this->myBrick->getNumExternalNodes();
  
  if ( numExtNodes != 8) {
    return false;
  }
  
  // Calculate signed volumes for point and four faces of brick
  // in order for the point to be included they must all be 
  // of negative sign given the orientations used
  // also global X+ is 1-4
  // global Y+ is  2-3
  // global Z+ is 1-5
  // implicit is the assupmtion that element is created
  // as such that lower left node is node 1, total = 8
  // bottom plane is nodes : 1-2-3-4
  // top plane is nodes : 5-6-7-8
  // and this is the order by which the element is 
  // defined


  // must be less than 0
  double bFSV = this->evalSignedVol(
				    this->myDomain->getNode(extNodes(0))->getCrds(),
				    this->myDomain->getNode(extNodes(1))->getCrds(),
				    this->myDomain->getNode(extNodes(2))->getCrds(),
				    pP);
  // must be less than 0
  double tFSV = this->evalSignedVol(
				    this->myDomain->getNode(extNodes(6))->getCrds(),
				    this->myDomain->getNode(extNodes(5))->getCrds(),
				    this->myDomain->getNode(extNodes(4))->getCrds(),
				    pP);
  // must be less than 0
  double lFSV = this->evalSignedVol(
				    this->myDomain->getNode(extNodes(3))->getCrds(),
				    this->myDomain->getNode(extNodes(7))->getCrds(),
				    this->myDomain->getNode(extNodes(4))->getCrds(),
				    pP);
  // must be less than 0 
  double rFSV = this->evalSignedVol(
				    this->myDomain->getNode(extNodes(2))->getCrds(),
				    this->myDomain->getNode(extNodes(1))->getCrds(),
				    this->myDomain->getNode(extNodes(5))->getCrds(),
				    pP);
  if ( (bFSV < 0) && (tFSV < 0) && (lFSV < 0) 
       && (rFSV < 0)) {
    return true;
  }

  return false;
}


double
GeometricBrickDecorator::evalSignedVol(const Vector &pA, 
				       const Vector &pB,
				       const Vector &pC,
				       const Vector &pP)
{
  return (1.0/6.0*(
	  - this->eval3OrderDet(pB, pC, pP) 
	  + this->eval3OrderDet(pA, pC, pP) 
	  - this->eval3OrderDet(pA, pB, pP)
	  + this->eval3OrderDet(pA, pB, pC)
	  ));
}


double
GeometricBrickDecorator::eval3OrderDet(const Vector &pA, 
				       const Vector &pB,
				       const Vector &pC)
{
  return (
	  pA(0)*pB(1)*pC(2) + 
	  pA(1)*pB(2)*pC(0) +
	  pA(2)*pB(0)*pC(1) -
	  pA(0)*pB(2)*pC(1) -
	  pA(1)*pB(0)*pC(2) -
	  pA(2)*pB(1)*pC(0) );
}

void
GeometricBrickDecorator::getFace(int which, ID &face, ID& faceID)
{
  const ID extNodes = this->myBrick->getExternalNodes();
  /*
    which == 1 -> bottom
    which == 2 -> top
    which == 3 -> left
    which == 4 -> right
    which == 5 -> front
    which == 6 -> rear
  */

  switch (which) {
  case 1:
    face(0) = extNodes(0);
    face(1) = extNodes(1);
    face(2) = extNodes(2);
    face(3) = extNodes(3);
    faceID(0) = 0;
    faceID(1) = 1;
    faceID(2) = 2;
    faceID(3) = 3;
    break;
  case 2:
    face(0) = extNodes(4);
    face(1) = extNodes(5);
    face(2) = extNodes(6);
    face(3) = extNodes(7);
    faceID(0) = 4;
    faceID(1) = 5;
    faceID(2) = 6;
    faceID(3) = 7;
    break;
  case 3:
    face(0) = extNodes(3);
    face(1) = extNodes(0);
    face(2) = extNodes(4);
    face(3) = extNodes(7);
    faceID(0) = 3;
    faceID(1) = 0;
    faceID(2) = 4;
    faceID(3) = 7;
    break;
  case 4:
    face(0) = extNodes(2);
    face(1) = extNodes(1);
    face(2) = extNodes(5);
    face(3) = extNodes(6);
    faceID(0) = 2;
    faceID(1) = 1;
    faceID(2) = 5;
    faceID(3) = 6;
    break;
  case 5:
    face(0) = extNodes(0);
    face(1) = extNodes(1);
    face(2) = extNodes(5);
    face(3) = extNodes(4);
    faceID(0) = 0;
    faceID(1) = 1;
    faceID(2) = 5;
    faceID(3) = 4;
    break;
  case 6:
    face(0) = extNodes(3);
    face(1) = extNodes(2);
    face(2) = extNodes(6);
    face(3) = extNodes(7);
    faceID(0) = 3;
    faceID(1) = 2;
    faceID(2) = 6;
    faceID(3) = 7;
    break;
  default:
    face(0) = -1;
    face(1) = -1;
    face(2) = -1;
    face(3) = -1;
    faceID(0) = -1;
    faceID(1) = -1;
    faceID(2) = -1;
    faceID(3) = -1;
    std::cerr << " ERROR in GeometricBrickDecorator L.233 \n";
    break;
  }
}

double
GeometricBrickDecorator::getMinMaxCrds(int which, int whichtoret) 
{
  /* 
     which goes from 1-> 3
  */

  double ret;

  Node** nodes = this->myBrick->getNodePtrs();
  double maxX = ((Vector&) nodes[1]->getCrds())(0);
  double minX = ((Vector&) nodes[0]->getCrds())(0);
  double maxY = ((Vector&) nodes[2]->getCrds())(1);
  double minY = ((Vector&) nodes[0]->getCrds())(1);
  double maxZ = ((Vector&) nodes[4]->getCrds())(2);
  double minZ = ((Vector&) nodes[0]->getCrds())(2);

  switch(which) {
  case 1:
    if (whichtoret<=0) {
      ret=minX;
    }
    else {
      ret = maxX;
    }
    break;
  case 2:
    if (whichtoret<=0) {
      ret=minY;
    }
    else {
      ret = maxY;
    }
    break;
  case 3:
    if (whichtoret<=0) {
      ret=minZ;
    }
    else {
      ret = maxZ;
    }
    break;
  default:
    std::cout << " ERROR ERROR ERROR in geometric brick decorator L.252 \n";
    ret = 0.0;
  }
  
  return ret;
}

bool
GeometricBrickDecorator::compareFaceToFace(int which, ID &faceOther) 
{
  ID faceID(4);
  ID myFace(4);
  this->getFace(which, myFace, faceID);
  return ( (myFace(0) == faceOther(0)) &&
	   (myFace(1) == faceOther(1)) &&
	   (myFace(2) == faceOther(2)) &&
	   (myFace(3) == faceOther(3)) 
	   );
}

bool
GeometricBrickDecorator::isFaceinPlane(int which, const Vector &pP) 
{
  /*
    general form to be implemented in the future 
  */
  return false;
}



bool
GeometricBrickDecorator::isFaceinVertPlane(int which, double xy, double zmin, double zmax, int whichCrd) 
{
  /*
    which -> which face : look at ::getFace(...)
    whichCrd -> 3 for Z
                2 for Y
		1 for X
  */
  
  int crd = 2;
  if ( (which == 1) || (which == 2) ) {
    crd = 3;
  }
  if ( (which == 3) || (which == 4) ) {
    crd = 1;
  }

  ID face(4);
  ID faceID(4);
  getFace(which, face, faceID);
  
  double maxE = getMinMaxCrds(whichCrd, 1);
  
  double minE = getMinMaxCrds(whichCrd, -1);
  Node** nodes = this->myBrick->getNodePtrs();
  Node* ndptr;
  double d0 = 0.0; 
  ndptr = nodes[faceID(0)];
  if (ndptr == 0 ) 
    opserr << " severe error NULL node ptr GeomDec L.294 \n" ;
  d0 = (ndptr->getCrds())(crd-1);
  double d1 = 0.0; 
  ndptr = nodes[faceID(1)];  
  if (ndptr == 0 ) 
    opserr << " severe error NULL node ptr GeomDec L.299 \n" ;
  d1 = (ndptr->getCrds())(crd-1);
  double d2 = 0.0; 
  ndptr = nodes[faceID(2)];  
  if (ndptr == 0 ) 
    opserr << " severe error NULL node ptr GeomDec L.304 \n" ;
  d2 = (ndptr->getCrds())(crd-1);
  double d3 = 0.0; 
  ndptr = nodes[faceID(3)];  
  if (ndptr == 0 ) 
    opserr << " severe error NULL node ptr GeomDec L.294 \n" ;
  d3 = (ndptr->getCrds())(crd-1);

  return (( d0== xy) && (d1 == xy) && ( d2 == xy) && (d3 == xy) &&
	  (maxE <= zmax) &&
	  (minE >= zmin) 
	  );
    
}

/* 
   Private methods used for the publc methods to identify the 
   DRM boundary layer elements
*/


bool 
GeometricBrickDecorator::isLeftBoundary(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
  return (
	  isFaceinVertPlane(3, xMin, zMin, zMax, 3) &&
	  isFaceinVertPlane(3, xMin, yMin, yMax, 2)
	  );
}

bool 
GeometricBrickDecorator::isRightBoundary(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
  return (
	  isFaceinVertPlane(4, xMax, zMin, zMax, 3) &&
	  isFaceinVertPlane(4, xMax, yMin, yMax, 2)
	  );
}

bool 
GeometricBrickDecorator::isFrontBoundary(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
  return (
	  isFaceinVertPlane(5, yMin, zMin, zMax, 3) &&
	  isFaceinVertPlane(5, yMin, xMin, xMax, 1)
	  );
}

bool 
GeometricBrickDecorator::isRearBoundary(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
  return (
	  isFaceinVertPlane(6, yMax, zMin, zMax, 3) &&
	  isFaceinVertPlane(6, yMax, xMin, xMax, 1)
	  );
}

bool 
GeometricBrickDecorator::isBottomBoundary(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
  return (
	  isFaceinVertPlane(1, zMin, xMin, xMax, 1) &&
	  isFaceinVertPlane(1, zMin, yMin, yMax, 2) 
	  );
}

bool
GeometricBrickDecorator::isBoundaryLayerEle(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
  if (	isLeftBoundary(xMin, xMax, yMin, yMax, zMin, zMax) ||

	isRightBoundary(xMin, xMax, yMin, yMax, zMin, zMax) ||
	
	isFrontBoundary(xMin, xMax, yMin, yMax, zMin, zMax) ||
	
	isRearBoundary(xMin, xMax, yMin, yMax, zMin, zMax) ||

	isBottomBoundary(xMin, xMax, yMin, yMax, zMin, zMax) 
	) 
    {
      return true;
    }
  return false;
}
