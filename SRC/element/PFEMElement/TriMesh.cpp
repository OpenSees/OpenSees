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
                                                                        
// $Revision: 1.0 $
// $Date: 2015-10-19 $
                                                                        
// Written: Minjie Zhu
// Created: Oct. 19
//
// Description: The class TriMesh is for meshing triangles
//

#include "TriMesh.h"
#include <Node.h>
#include <NodeIter.h>
#include <Domain.h>
#include <cmath>
#include <MeshRegion.h>
#include <string>
#include "TriangleMeshGenerator.h"
#include <elementAPI.h>
#include <Element.h>

extern int OPS_PFEMElement2D(Domain& theDomain, const ID& elenodes, ID& eletags);
extern int OPS_PFEMElement2DCompressible(Domain& theDomain, const ID& elenodes, ID& eletags);
extern int OPS_PFEMElement2DBubble(Domain& theDomain, const ID& elenodes, ID& eletags);
extern int OPS_Tri31(Domain& theDomain, const ID& elenodes, ID& eletags);

int OPS_TriMesh(Domain& domain)
{
    if (OPS_GetNumRemainingInputArgs() < 4) {
	opserr<<"WARNING: want rtag? ndf? size? numnodes? nd1? ... bound1? ... eletype? ...\n";
	return -1;
    }

    int num = 2;
    int rtag[2];
    if (OPS_GetIntInput(&num,&rtag[0]) < 0) {
	opserr<<"WARNING: failed to read region tag\n";
	return -1;
    }

    num = 1;
    double size;
    if (OPS_GetDoubleInput(&num,&size) < 0) {
	opserr<<"WARNING: failed to read mesh size\n";
	return -1;
    }

    num = 1;
    int numnodes;
    if (OPS_GetIntInput(&num,&numnodes) < 0) {
	opserr<<"WARNING: failed to read number of nodes\n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() < 2*numnodes) {
	opserr<<"WARNING: insufficient input args\n";
	return -1;
    }

    ID nodes(numnodes);
    if (OPS_GetIntInput(&numnodes, &nodes(0)) < 0) {
	opserr<<"WARNING: failed to read node tags\n";
	return -1;
    }

    ID bound(numnodes);
    if (OPS_GetIntInput(&numnodes, &bound(0)) < 0) {
	opserr<<"WARNING: failed to read boundary values\n";
	return -1;
    }

    TriMesh mesh(domain,rtag[1]);

    return mesh.mesh(rtag[0],size,nodes,bound);
}

int OPS_TriReMesh(Domain& domain, int ndf)
{
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr<<"WARNING: want alpha? numfreeregions? rtag1? ... numfixregions? rtag1? ...\n";
	return -1;
    }

    int num = 1;
    double alpha;
    if (OPS_GetDoubleInput(&num,&alpha) < 0) {
	opserr<<"WARNING: failed to read alpha\n";
	return -1;
    }

    int numfreeregions;
    if (OPS_GetIntInput(&num,&numfreeregions) < 0) {
	opserr<<"WARNING: failed to read number of regions to be remeshed\n";
	return -1;
    }
    if (numfreeregions < 1) {
	opserr << "WARNING: numfreeregions < 1\n";
	return -1;
    }
    
    if (OPS_GetNumRemainingInputArgs() < numfreeregions+1) {
	opserr<<"WARNING: insufficient args for free regions\n";
	return -1;
    }

    ID rtagsfree(numfreeregions);
    if (OPS_GetIntInput(&numfreeregions, &rtagsfree(0)) < 0) {
	opserr<<"WARNING: failed to read region tags\n";
	return -1;
    }

    ID rtagsfix;
    if (OPS_GetNumRemainingInputArgs() > 0) {
	int numfixregions;
	if (OPS_GetIntInput(&num,&numfixregions) < 0) {
	    opserr<<"WARNING: failed to read number of regions to be fixed\n";
	    return -1;
	}
	
	if (numfixregions > 0) {
	    if (OPS_GetNumRemainingInputArgs() < numfixregions) {
		opserr<<"WARNING: insufficient args for fixed regions\n";
		return -1;
	    }
	    
	    rtagsfix.resize(numfixregions);
	    if (OPS_GetIntInput(&numfixregions, &rtagsfix(0)) < 0) {
		opserr<<"WARNING: failed to read region tags\n";
		return -1;
	    }
	}
    }


    // element args
    if (OPS_GetNumRemainingInputArgs() == 0) {
	opserr<<"WARNING: element args are not given\n";
	return -1;
    }

    TriMesh mesh(domain,ndf);

    return mesh.mesh(alpha,rtagsfree,rtagsfix);
}



TriMesh::TriMesh(Domain& domain, int n)
    :theDomain(&domain), ndf(n)
{
}

TriMesh::~TriMesh()
{
}

int
TriMesh::mesh(int rtag, double size, const ID& nodes,const ID& bound)
{
    // check
    if(size <= 0) {
	opserr<<"WARNING: mesh size <= 0\n";
	return -1;
    }
    if(nodes.Size() < 3) {
	opserr<<"WARNING: input number of nodes < 3\n";
	return -1;
    }
    if(bound.Size() < nodes.Size()) {
	opserr<<"WARNING: the number of boundary ID < number of nodes\n";
	return -1;
    }
    
    // calling mesh generator
    TriangleMeshGenerator gen;
    for(int i=0; i<nodes.Size(); i++) {
	// get node
	Node* theNode = theDomain->getNode(nodes(i));
	if(theNode == 0) {
	    opserr<<"WARNING: node "<<nodes(i)<<" is not defined\n";
	    return -1;
	}
	const Vector& crds = theNode->getCrds();
	if(crds.Size() < 2) {
	    opserr<<"WARNING: ndm < 2\n";
	    return -1;
	}
	// add point
	gen.addPoint(crds(0), crds(1));

	// add segment
	int p1 = i;
	int p2;
	if(i==nodes.Size()-1) {
	    p2 = 0;
	} else {
	    p2 = i+1;
	}
	if(bound(i) == 0) {
	    gen.addSegment(p1,p2,-1);
	} else {
	    gen.addSegment(p1,p2,0);
	}
    }

    // meshing
    gen.mesh(size*size*0.5);

    // get node tag
    NodeIter& theNodes = theDomain->getNodes();
    Node* theNode = theNodes();
    int currtag = 0;
    if(theNode != 0) currtag = theNode->getTag();

    // get nodes
    ID regnodes(0,gen.getNumPoints());
    ID ptmark(gen.getNumPoints());
    ID ptnode(gen.getNumPoints());
    int index = 0;
    for(int i=0; i<nodes.Size(); i++) {
	int j = i-1;
	if(i==0) j = nodes.Size()-1;
	if(bound(i)!=0 && bound(j)!=0) {
	    regnodes[index++] = nodes(i);
	    ptmark(i) = 0;
	} else {
	    ptmark(i) = -1;
	}
	ptnode(i) = nodes(i);
    }
    for(int i=nodes.Size(); i<gen.getNumPoints(); i++) {
	// get point
	double x, y;
	int mark = 0;
	gen.getPoint(i,x,y,mark);

	// if on unwanted boundary
	if(mark == -1) {
	    ptmark(i) = -1;
	    continue;
	} else {
	    ptmark(i) = 0;
	}

	// create node
	Node* node = new Node(--currtag,ndf,x,y);
	if(node == 0) {
	    opserr<<"run out of memory for creating Node\n";
	    return -1;
	}
	if(theDomain->addNode(node) == false) {
	    opserr<<"Failed to add node to domain\n";
	    delete node;
	    return -1;
	}

	// add to region
	regnodes[index++] = currtag;
	ptnode(i) = currtag;
    }

    // get elenodes
    ID regelenodes(0,gen.getNumTriangles());
    index = 0;
    for(int i=0; i<gen.getNumTriangles(); i++) {
	int p1,p2,p3;
	gen.getTriangle(i,p1,p2,p3);
	if(ptmark(p1)==0 && ptmark(p2)==0 && ptmark(p3)==0) {
	    regelenodes[index++] = ptnode(p1);
	    regelenodes[index++] = ptnode(p2);
	    regelenodes[index++] = ptnode(p3);
	}
    }

    // get region
    MeshRegion* theRegion = theDomain->getRegion(rtag);
    if(theRegion == 0) {
	theRegion = new MeshRegion(rtag);
	if(theDomain->addRegion(*theRegion) < 0) {
	    opserr<<"WARNING: failed to add region\n";
	    return -1;
	}
    }

    // add to region
    theRegion->setNodes(regnodes);
    
    return 0;
}

int
TriMesh::mesh(double alpha, const ID& rtags2, const ID& rtags)
{
    // get all regions
    int numregions = rtags.Size()+rtags2.Size();
    int numpoints = 0;
    std::vector<MeshRegion*> regions(numregions);
    for(int i=0; i<rtags.Size(); i++) {
	MeshRegion* region = theDomain->getRegion(rtags(i));
	if(region == 0) {
	    opserr<<"WARNING: region "<<rtags(i)<<" is not defined\n";
	    return -1;
	}
	numpoints += region->getNodes().Size();
	regions[i] = region;
    }
    for(int i=0; i<rtags2.Size(); i++) {
	MeshRegion* region = theDomain->getRegion(rtags2(i));
	if(region == 0) {
	    opserr<<"WARNING: region "<<rtags2(i)<<" is not defined\n";
	    return -1;
	}
	numpoints += region->getNodes().Size();
	regions[i+rtags.Size()] = region;

	// remove elements
	const ID& eles = region->getExtraEles();
	for (int j=0; j<eles.Size(); j++) {
	    Element* ele = theDomain->removeElement(eles(j));
	    if (ele != 0) {
		delete ele;
	    }
	}

    }

    // get all points
    ID ptreg(0,numpoints), ptnode(0,numpoints);
    for(int i=0; i<numregions; i++) {
	const ID& nodes = regions[i]->getNodes();
	for(int j=0; j<nodes.Size(); j++) {
	    int index = ptreg.Size();
	    ptreg[index] = i;
	    ptnode[index] = nodes(j);
	}
    }

    // calling mesh generator
    TriangleMeshGenerator gen;
    for(int i=0; i<ptnode.Size(); i++) {
	// get node
	Node* theNode = theDomain->getNode(ptnode(i));
	if(theNode == 0) {
	    opserr<<"WARNING: node "<<ptnode(i)<<" is not defined\n";
	    return -1;
	}
	const Vector& crds = theNode->getCrds();
	const Vector& disp = theNode->getTrialDisp();
	if(crds.Size() < 2 || disp.Size() < 2) {
	    opserr<<"WARNING: ndm < 2 or ndf < 2\n";
	    return -1;
	}
	
	// add point
	gen.addPoint(crds(0)+disp(0), crds(1)+disp(1));
    }

    // meshing
    gen.remesh(alpha);

    // get elenodes
    std::vector<ID> regelenodes(rtags2.Size());
    for(int i=0; i<gen.getNumTriangles(); i++) {
	int p1,p2,p3;
	gen.getTriangle(i,p1,p2,p3);

	// if all connected to fixed regions
	int numfix = rtags.Size();
	if(ptreg(p1)<numfix && ptreg(p2)<numfix && ptreg(p3)<numfix) {
	    continue;
	}
	
	// which region it belongs to
	int reg = ptreg(p1);
	if(reg<numfix || (reg>ptreg(p2) && ptreg(p2)>=numfix)) reg = ptreg(p2);
	if(reg<numfix || (reg>ptreg(p3) && ptreg(p3)>=numfix)) reg = ptreg(p3);
	reg -= numfix;

	// elenodes
	int index = regelenodes[reg].Size();
	regelenodes[reg][index++] = ptnode(p1);
	regelenodes[reg][index++] = ptnode(p2);
	regelenodes[reg][index++] = ptnode(p3);
    }

    // all elenodes
    ID allelenodes;
    for(int i=0; i<(int)regelenodes.size(); i++) {
	int num = allelenodes.Size();
	int num1 = regelenodes[i].Size();
	allelenodes.resize(num+num1);
	for (int j=0; j<num1; j++) {
	    allelenodes(num+j) = regelenodes[i](j);
	}
    }

    // create elements
    std::string eletype = OPS_GetString();
    ID eletags;
    if (eletype=="PFEMElement2D") {
	if (OPS_PFEMElement2D(*theDomain,allelenodes,eletags) < 0) {
	    return -1;
	}
    } else if (eletype=="PFEMElement2DQuasi") {
	if (OPS_PFEMElement2DCompressible(*theDomain,allelenodes,eletags) < 0) {
	    return -1;
	}
    } else if (eletype=="PFEMElement2DBubble") {
	if (OPS_PFEMElement2DBubble(*theDomain,allelenodes,eletags) < 0) {
	    return -1;
	}

    }  else if (eletype=="tri31") {
	if (OPS_Tri31(*theDomain,allelenodes,eletags) < 0) {
	    return -1;
	}

    } else {
	opserr<<"WARNING: element "<<eletype.c_str()<<" can't be used for tri mesh at this time\n";
	return -1;
    }

    // add to region
    int num = 0;

    for(int i=rtags.Size(); i<numregions; i++) {
	int num1 = regelenodes[i-rtags.Size()].Size()/3;
	ID eletag(num1);
	for (int j=0; j<num1; j++) {
	    eletag(j) = eletags(num+j);
	}
	num += num1;
	regions[i]->setExtraEles(eletag);
    }
    return 0;
}
