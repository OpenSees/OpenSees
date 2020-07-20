#include "Mesh3DSubdomain.h"
#include <stdlib.h>
#include <iostream>
#include <OPS_Globals.h>
#include <StandardStream.h>
#include <Timer.h>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <NDMaterial.h>

#include <ArrayOfTaggedObjects.h>

// includes for the domain classes
#include <Domain.h>
#include <ShadowSubdomain.h>
#include <ActorSubdomain.h>
#include <Node.h>
#include <Element.h>
#include <ZeroLength.h>
#include <ElasticMaterial.h>
#include <ElasticIsotropicThreeDimensional.h>

#include <SP_Constraint.h>

#include <SingleDomAllSP_Iter.h>
#include <NodeIter.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <Graph.h>
#include <Vector.h>
#include <Matrix.h>
#include "GeometricBrickDecorator.h"
#include <math.h>
using namespace std;


#include <Element.h>
#include <Node.h>
#include <DummyNode.h>
#include <RigidBeam.h>
#include <ElementIter.h>
#include <LinearCrdTransf3d.h>
#include <ElasticBeam3d.h>
#include <MP_Constraint.h>
#include <classTags.h>
//general constructor
Mesh3DSubdomain::Mesh3DSubdomain(Domain * inpDomain) 
{

  this->myDomain = inpDomain;
  this->myLastEle = 0;
  this->myLastNode = 0;
  this->myLastStartNodeTag = 0;
}

//general constructor
Mesh3DSubdomain::Mesh3DSubdomain() 
{

  this->myDomain = NULL;
  this->myLastEle = 0;
  this->myLastNode = 0;
  this->myLastStartNodeTag = 0;
}

//destructor
Mesh3DSubdomain::~Mesh3DSubdomain()
{ 
}

int
Mesh3DSubdomain::lastEle()
{
  return this->myLastEle;
}

int
Mesh3DSubdomain::lastNode()
{
  return this->myLastNode;
}

int
Mesh3DSubdomain::lastStartNode()
{
  return this->myLastStartNodeTag;
}



/*
  Currently works for subdomains whose top surface is a square so that
  it is possible to limit the application of the lysmer elements to all
  the exterior nodes but the ones at the top surface
*/

void 
Mesh3DSubdomain::allocate_e_Nodes(double xMin, double xMax, 
				  double yMin, double yMax,
				  double zMin, double zMax,
				  std::map<int,int>& eNodes)
{
  Vertex *theVertex;
  double xCoord, yCoord, zCoord;
  Graph  & nodeGraph = this->myDomain->getNodeGraph();
  VertexIter  & theVertices = nodeGraph.getVertices();
  int count = 0;
  bool enode_found = false;
  while ((theVertex = theVertices() ) != 0) {
    int tag = theVertex->getRef();

    const Vector & coords = this->myDomain->getNode(tag)->getCrds();
    
    xCoord = coords(0);
    yCoord = coords(1);
    zCoord = coords(2);
    
    if ( ((xCoord == xMin) || (xCoord == xMax)) &&
	 (yCoord >= yMin) && (yCoord<=yMax) &&
	 (zCoord >= zMin) && (zCoord<=zMax)) 
      {
	enode_found = true;
      }

    if ( ((yCoord == yMin) || (yCoord == yMax)) &&
	 (xCoord >= xMin) && (xCoord<=xMax) &&
	 (zCoord >= zMin) && (zCoord<=zMax)) 
      {
	enode_found = true;
      }

    if ( ((zCoord == zMin)) && 
	 (((xCoord >= xMin) && (xCoord<=xMax)) &&
	  ((yCoord >= yMin) && (yCoord<=yMax)) ))
      {
	enode_found = true;
      }
    if ((zCoord == zMax) &&
	((xCoord == xMin) || (xCoord ==xMax)) &&
	((yCoord >= yMin) && (yCoord <= yMax)))
      { 
	enode_found = true;
      }
    if ((zCoord == zMax) &&
	((yCoord == yMin) || (yCoord==yMax)) &&
	((xCoord >= xMin) && (xCoord <= xMax)))
      {
	enode_found = true;
      }
    if (enode_found) {
      eNodes[tag] = count;
      count++;
      enode_found = false;
    }
  }
}

void
Mesh3DSubdomain::allocateBoundaryLayerElements(double xMin, double xMax,
					       double yMin, double yMax,
					       double zMin, double zMax,
					       std::map<int,Element*>& elements,
					       std::map<int,Vector*>& storage,
					       std::map<int,int>& storage2)
{
  Vertex* theVertex;
//  int count =0;
//  bool ele_found = false;
  GeometricBrickDecorator* myHelper = new GeometricBrickDecorator();
  myHelper->setDomain(this->myDomain);
  Graph& eleGraph = this->myDomain->getElementGraph();
  VertexIter& theVertices = eleGraph.getVertices();
  while ((theVertex = theVertices() ) != 0) {
    int tag = theVertex->getRef();
    Element* ele = this->myDomain->getElement(tag);
    if ((ele->getClassTag() == ELE_TAG_EightNode_Brick_u_p) || (ele->getClassTag() == ELE_TAG_Brick) || (ele->getClassTag() == ELE_TAG_FLBrick) ) {
      myHelper->setBrick(ele);
      if ( myHelper->isBoundaryLayerEle(xMin, xMax, yMin, yMax, zMin, zMax) ) {
	Vector* ptr = new Vector(24);
	ptr->Zero();
	storage[tag] = ptr; 
	storage2[tag] = -1;
	elements[tag] = ele;
      }
    }
  }
  delete myHelper;
}
