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
// Description: The class LineMesh is for meshing lines.
//

#include "LineMesh.h"
#include <Node.h>
#include <NodeIter.h>
#include <Domain.h>
#include <cmath>
#include <MeshRegion.h>
#include <string>
#include <elementAPI.h>

extern int OPS_ElasticBeam2d(Domain& theDomain, const ID& elenodes, ID& eletags);
extern int OPS_ForceBeamColumn2d(Domain& theDomain, const ID& elenodes, ID& eletags);
extern int OPS_DispBeamColumn2d(Domain& theDomain, const ID& elenodes, ID& eletags);

int OPS_LineMesh(Domain& domain, int ndm)
{
    if (OPS_GetNumRemainingInputArgs() < 4) {
	opserr<<"WARNING: want rtag? ndf? size? numnodes? nd1? ... bound1? ... eletype? ...\n";
	return -1;
    }

    // get inputs
    int num = 2;
    int rtag[2];
    if (OPS_GetIntInput(&num,rtag) < 0) {
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

    // create mesher
    LineMesh mesh(domain,ndm,rtag[1]);

    // mesh
    if (mesh.mesh(rtag[0],size,nodes,bound) < 0) {
	opserr<<"WARNING: failed to mesh line\n";
	return -1;
    }

    return 0;
}

LineMesh::LineMesh(Domain& domain, int m, int n)
    :theDomain(&domain), ndm(m), ndf(n)
{
}

LineMesh::~LineMesh()
{
}


int
LineMesh::mesh(int rtag, double size, const ID& nodes,const ID& bound)
{
    if(size <= 0) {
	opserr<<"WARNING: mesh size <= 0\n";
	return -1;
    }
    
    int nodesize = nodes.Size();
    if(nodesize < 2) {
	opserr<<"WARNING: the number of input nodes < 2\n";
	return -1;
    }
    if(bound.Size() < nodesize) {
	opserr<<"WARNING: the number of boundary ID < number of nodes\n";
	return -1;
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
    
    // get region nodes
    ID regnodes, regelenodes;
    for(int i=0; i<nodesize; i++) {

	// last node
	int num = regnodes.Size();
	if (i == nodesize-1) {
	    if (bound(i) != 0) {
		regnodes[num] = nodes(i);
	    }
	    break;
	}

	// end nodes
	int nd1 = nodes(i);
	int nd2 = nodes(i+1);

	// mesh line
	ID lnodes, elenodes;
	if(this->mesh(size,nd1,nd2,lnodes,elenodes) < 0) {
	    opserr<<"WARNING: failed to mesh line\n";
	    return -1;
	}
	
	// add nodes to region
	if(bound(i) != 0) {
	    regnodes.resize(num+lnodes.Size()+1);
	    regnodes(num) = nd1;
	    num++;
	} else {
	    regnodes.resize(num+lnodes.Size());
	}
	for(int j=0; j<lnodes.Size(); j++) {
	    regnodes(num+j) = lnodes(j);
	}
	
	// add elenodes to region
	num = regelenodes.Size();
	regelenodes.resize(num+elenodes.Size());
	for(int j=0; j<elenodes.Size(); j++) {
	    regelenodes(num+j) = elenodes(j);
	}
    }

    // add to region
    theRegion->setNodes(regnodes);
        
    // create elements
    if (OPS_GetNumRemainingInputArgs() > 0) {
	
	std::string eletype = OPS_GetString();
	ID eletags;
	if (eletype=="elasticBeamColumn") {
	    if (ndm == 2) {
		if (OPS_ElasticBeam2d(*theDomain,regelenodes,eletags) < 0) {
		    return -1;
		}
	    }
	} else if (eletype == "forceBeamColumn") {
	    if (ndm == 2) {
		if (OPS_ForceBeamColumn2d(*theDomain,regelenodes,eletags) < 0) {
		    return -1;
		}
	    }

	} else if (eletype == "dispBeamColumn") {
	    if (ndm == 2) {
		if (OPS_DispBeamColumn2d(*theDomain,regelenodes,eletags) < 0) {
		    return -1;
		}
	    }

	} else {
	    opserr<<"WARNING: element "<<eletype.c_str()<<" can't be used for line mesh at this time\n";
	    return -1;
	}
	theRegion->setExtraEles(eletags);
    }

    

    return 0;
}

int
LineMesh::mesh(double size, int tag1, int tag2, ID& nodes, ID& elenodes)
{
    if(tag1 == tag2) {
	opserr<<"WARNING: same tags are given as ends of line\n";
	return -1;
    }
    
    // get end nodes
    Node* nd1 = theDomain->getNode(tag1);
    Node* nd2 = theDomain->getNode(tag2);
    if(nd1==0) {
	opserr<<"node "<<tag1<<" is not definde\n";
	return -1;
    }
    if(nd2==0) {
	opserr<<"node "<<tag2<<" is not definde\n";
	return -1;
    }

    // node coordinats
    const Vector& crds1 = nd1->getCrds();
    const Vector& crds2 = nd2->getCrds();
    if(crds1.Size()!=ndm || crds2.Size()!=ndm) {
	opserr<<"WARNING: node "<<tag1<<" and "<<tag2;
	opserr<<" have different ndm to the model\n";
	return -1;
    }

    // num of ele
    Vector incr = crds2-crds1;
    int nele = ceil(incr.Norm()/size);
    if(nele == 0) {
	opserr<<"WARNING: two nodes are at same location\n";
	return -1;
    }

    // increment
    incr /= nele;

    // get node tag
    NodeIter& theNodes = theDomain->getNodes();
    Node* theNode = theNodes();
    int currtag = 0;
    if(theNode != 0) currtag = theNode->getTag();

    // creat nodes
    if(nele > 1) {
	nodes.resize(nele-1);
    }
    Vector crds = crds1;
    for(int i=0; i<nele-1; i++) {
	crds += incr;
	Node* node = 0;
	if(ndm == 1) {
	    node = new Node(--currtag,ndf,crds(0));
	} else if(ndm == 2) {
	    node = new Node(--currtag,ndf,crds(0),crds(1));
	} else if(ndm == 3) {
	    node = new Node(--currtag,ndf,crds(0),crds(1),crds(2));
	}
	if(node == 0) {
	    opserr<<"run out of memory for creating Node\n";
	    return -1;
	}
	if(theDomain->addNode(node) == false) {
	    opserr<<"Failed to add node to domain\n";
	    delete node;
	    return -1;
	}
	nodes(i) = currtag;
    }

    // line elements
    int index = 0;
    elenodes.resize(nele*2);
    elenodes(index++) = tag1;
    if(nele == 1) {
	elenodes(index++) = tag2;
    } else {
	elenodes(index++) = nodes(0);
    }
    for(int i=1; i<nele-1; i++) {
	elenodes(index++) = nodes(i-1);
	elenodes(index++) = nodes(i);
    }
    if(nele > 1) {
	elenodes(index++) = nodes(nele-2);
	elenodes(index++) = tag2;
    }

    return 0;
}
