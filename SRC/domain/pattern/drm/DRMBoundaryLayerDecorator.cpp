#include "DRMBoundaryLayerDecorator.h"
#include <Brick.h>

DRMBoundaryLayerDecorator::
DRMBoundaryLayerDecorator()
{
  std::map<int,int> mapa;
  std::set<int,std::less<int > > seta;

  this->myDomain = NULL;
  this->myBrick = NULL;
  this->eNodeArray = NULL;
  this->eNodeID = NULL;
  this->eNodeMap = mapa;
  this->eNodeSet = seta;
  this->cons = -1;

  Peff_k = new Vector(24);
  Peff_d = new Vector(24);
  Peff_m = new Vector(24);
}


DRMBoundaryLayerDecorator::
~DRMBoundaryLayerDecorator()
{
  delete Peff_k;
  delete Peff_d;
  delete Peff_m;
}

void 
DRMBoundaryLayerDecorator::setArray(int *eNodes)
{
  this->eNodeArray = eNodes;
  this->cons = 0;
}

void 
DRMBoundaryLayerDecorator::setID(ID* eNodes)
{
  this->eNodeID = eNodes;
  this->cons = 1;
}

void 
DRMBoundaryLayerDecorator::setMap(std::map<int,int> &eNodes)
{
  this->eNodeMap = eNodes;
  this->cons = 2;
}

void 
DRMBoundaryLayerDecorator::setSet(std::set<int,std::less<int > > &eNodes)
{
  this->eNodeSet = eNodes;
  this-> cons = 3;
}

void 
DRMBoundaryLayerDecorator::setDomain(Domain *theDomain)
{
  this->myDomain = theDomain;
}

void 
DRMBoundaryLayerDecorator::clearDomain()
{
  this->myDomain = NULL;
}

void 
DRMBoundaryLayerDecorator::setBrick(Element* brickTag)
{
  this->myBrick = brickTag;
}
  
void 
DRMBoundaryLayerDecorator::clearBrick()
{
  this->myBrick = NULL;
}

void 
DRMBoundaryLayerDecorator::clear()
{
  this->myDomain = NULL;
  this->myBrick = NULL;
}



void 
DRMBoundaryLayerDecorator::get_E_B_Nodes(ID &e, ID &b)
{
  /// Map version works for map having map[tag] = -1;

  const ID & nodes  = this->myBrick->getExternalNodes();
  int nodeTag;
  switch(this->cons) {
  case 0:
    for (int i=0; i<8; i++) {
      nodeTag = nodes(i);
      if (eNodeArray[nodeTag] == -1) {
	e(i) = -1;
      }
      else {
	b(i) = -1;
      }
    }
    break;
  case 1:
    {
      ID& tempNode = *this->eNodeID;
      for (int i=0; i<8; i++) {
	nodeTag = nodes(i);
	if (tempNode.getLocation(nodeTag) >= 0) {
	  e(i) = -1;
	}
	else {
	  b(i) = -1;
	}
      }
      break;
    }
  case 2:
    for (int i=0; i<8; i++) {
      nodeTag = nodes(i);
      if (!(eNodeMap.find(nodeTag) == eNodeMap.end() )) {
	e(i) = -1;
      }
      else {
	b(i) = -1;
      }
    }
    break;
  case 3:
    for (int i=0; i<8; i++) {
      nodeTag = nodes(i);
      if (eNodeSet.find(nodeTag) != eNodeSet.end()) {
	e(i) = -1;
      }
      else {
	b(i) = -1;
      }
    }
    break;
  default:
    break;
  }
}



void 
DRMBoundaryLayerDecorator::zeroSubmatrix(Matrix &M, int e, int b)
{
  for ( int row = 3*e; row<(3*e+3); row++) {
    for ( int col = 3*b; col<(3*b+3); col++) {
      M(row,col) = 0.0;
    }
  }
}

void 
DRMBoundaryLayerDecorator::computeDRMLoad(Vector &drmLoad, 
				      const Vector &displ,
				      const Vector &veloc,
				      const Vector &accel)
{
  Matrix *Ke = new Matrix(this->myBrick->getTangentStiff());
  if (Ke == 0) {
    opserr << " NO MATRIX Ke ALLOCATED \n";
  }
  Matrix *Ce = new Matrix(this->myBrick->getDamp());
  if (Ce == 0) {
     opserr << " NO MATRIX Ce ALLOCATED \n";
  }

  Matrix *Me = new Matrix(this->myBrick->getMass());
  if (Me == 0) {
     opserr << " NO MATRIX Me ALLOCATED \n";
  }

  ID e(8), b(8);
  e.Zero();
  b.Zero();
  
  this->get_E_B_Nodes(e, b);

  for (int i=0; i<8; i++) {
    for (int j=0; j<8; j++) {
      if (e(i) != b(j)) {
	zeroSubmatrix(*Ke, i, j);
	zeroSubmatrix(*Ce, i, j);
	zeroSubmatrix(*Me, i, j);
      }
    }
  }

  Peff_k->addMatrixVector(0.0, *Ke, displ, 1.0);
  Peff_d->addMatrixVector(0.0, *Ce, veloc, 1.0);
  Peff_m->addMatrixVector(0.0, *Me, accel, 1.0);
  
  for (int i=0; i<24; i++) {
    if (e(i/3) == -1) {
      drmLoad(i) = (*Peff_k)(i) + (*Peff_d)(i) + (*Peff_m)(i);
    }
    else {
      drmLoad(i) = -(*Peff_k)(i) - (*Peff_d)(i) -(*Peff_m)(i);
    }
  }

  delete Ke;
  delete Ce;
  delete Me;

}


void 
DRMBoundaryLayerDecorator::applyDRMLoad(double cfact,
					Vector &drmLoad,
					const Vector &displ,
					const Vector &veloc,
					const Vector &accel)
{
  Node *theNode;
  drmLoad.Zero();
  this->computeDRMLoad(drmLoad, displ, veloc, accel);
  Node** nodes  = this->myBrick->getNodePtrs();
  Vector load(3);
  for (int i=0; i<8; i++) {
    theNode = nodes[i];
    load.Zero();
    load(0) = cfact*drmLoad(i*3);
    load(1) = cfact*drmLoad(i*3+1);
    load(2) = cfact*drmLoad(i*3+2);
    theNode->addUnbalancedLoad(load);
  }
}
