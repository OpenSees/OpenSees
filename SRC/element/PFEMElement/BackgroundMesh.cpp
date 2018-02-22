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

// $Revision$
// $Date$
// $URL$
                                                                        
// Written: Minjie Zhu
//
// Description: This class defines the BackgroundMesh
//

#include <BackgroundMesh.h>
#include <Vector.h>
#include <Domain.h>
#include <elementAPI.h>
#include <ID.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <Matrix.h>
#include <NodeIter.h>
#include <cmath>
#include <Domain.h>
#include <Pressure_Constraint.h>
#include <Pressure_ConstraintIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <PFEMElement2DBubble.h>
#include <fstream>
#include <iostream>
#include <string.h>

void* OPS_PVDRecorder();

// static BackgroundMesh* backgroundMesh = new BackgroundMesh;
BackgroundMesh backgroundMesh;

double BackgroundMesh::pi = 3.141592653589793;

BackgroundMesh& OPS_GetBackgroundMesh()
{
    return backgroundMesh;
}

int OPS_BackgroundMesh()
{
    int ndm = OPS_GetNDM();

    BackgroundMesh* background = &OPS_GetBackgroundMesh();
    
    // check input
    if(OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"WARNING: no type is specified -- background\n";
	return -1;
    }
    if (ndm <1 || ndm > 3) {
	opserr<<"WARNING: bad ndm = "<<ndm<<"\n";
	return -1;
    }

    // type
    const char* type = OPS_GetString();

    // sub-type
    ID newnodes, elenodes, trinodes;
    if (strcmp(type,"meshsize") == 0) {

	if (OPS_GetNumRemainingInputArgs() <1) {
	    opserr<<"WARNING: insufficient args\n";
	    return -1;
	}

	// num divitions
	double size;
	int num = 1;
	if (OPS_GetDoubleInput(&num, &size) < 0) {
	    opserr<<"WARNING: failed to get mesh size\n";
	    return -1;
	}
	background->setMeshsize(size);

    } else if (strcmp(type, "pvd") == 0) {

	PVDRecorder* recorder = (PVDRecorder*)OPS_PVDRecorder();
	if (recorder == 0) {
	    opserr << "WARNING failed to create pvd recorder\n";
	    return -1;
	}
	background->setRecorder(recorder);

    } else if (strcmp(type, "record") == 0) {

	int ctag = 0;
	Domain* theDomain = OPS_GetDomain();
	if (theDomain == 0) return -1;

	double timestamp = theDomain->getCurrentTime();
	if (background->record(ctag, timestamp) < 0) {
	    opserr << "WARNING failed to record\n";
	    return -1;
	}

    } else if (strcmp(type,"mapToBack") == 0) {

	if (background->mapParticleToBack() < 0) return -1;

    } else if (strcmp(type, "mapFromBack") == 0) {

	if (background->mapBackToParticle() < 0) return -1;

    } else if (strcmp(type, "fix") == 0) {

	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr<<"WARNING: insufficient args\n";
	    return -1;
	}

	const char* dir = OPS_GetString();

	Vector min(ndm), max(ndm);

	double halfsize = background->getMeshsize()/2.0;

	if (strcmp(dir,"p") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < ndm) {
		opserr<<"WARNING: insufficient args\n";
		return -1;
	    }
	    if (OPS_GetDoubleInput(&ndm, &min(0)) < 0) {
		opserr<<"WARNING: failed to get lower left point\n";
		return -1;
	    }
	    for(int i=0; i<ndm; i++) {
		max(i) = min(i)+halfsize;
		min(i) -= halfsize;
	    }
	    
	} else if (strcmp(dir,"x") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 3) {
		opserr<<"WARNING: insufficient args\n";
		return -1;
	    }
	    Vector val(3);
	    int num = 3;
	    if (OPS_GetDoubleInput(&num, &val(0)) < 0) {
		opserr<<"WARNING: failed to get x position\n";
		return -1;
	    }
	    min(0) = val(0)-halfsize;
	    min(1) = val(1)-halfsize;
	    max(0) = val(0)+halfsize;
	    max(1) = val(2)+halfsize;
	    
	    
	} else if (strcmp(dir,"y") == 0) {

	    if (OPS_GetNumRemainingInputArgs() < 3) {
		opserr<<"WARNING: insufficient args\n";
		return -1;
	    }
	    Vector val(3);
	    int num = 3;
	    if (OPS_GetDoubleInput(&num, &val(0)) < 0) {
		opserr<<"WARNING: failed to get y position\n";
		return -1;
	    }
	    min(0) = val(1)-halfsize;
	    min(1) = val(0)-halfsize;
	    max(0) = val(2)+halfsize;
	    max(1) = val(0)+halfsize;
	} else {
	    opserr<<"WARNING: unidentified fix type\n";
	    return -1;
	}

	int num = OPS_GetNumRemainingInputArgs();
	ID fix(num);
	if (OPS_GetIntInput(&num, &fix(0)) < 0) {
	    opserr<<"WARNING: failed to get fix values\n";
	    return -1;
	}

	background->addFixInfo(min,max,fix);


    } else if (strcmp(type,"particle") == 0) {

	if (OPS_GetNumRemainingInputArgs() < 8) {
	    opserr<<"WARNING: insufficient args\n";
	    return -1;
	}

	// a geometry
	const char* geotype = OPS_GetString();
	Vector p1(ndm), p2(ndm), p3(ndm), p4(ndm);
	ID nump(2);
	if (strcmp(geotype,"quad") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 4*ndm+2) {
		opserr<<"WARNING: insufficient args\n";
		return -1;
	    }

	    // node coord
	    if (OPS_GetDoubleInput(&ndm, &p1(0)) < 0) {
		opserr<<"WARNING: failed to get cooridnates for first point\n";
		return -1;
	    }
	    if (OPS_GetDoubleInput(&ndm, &p2(0)) < 0) {
		opserr<<"WARNING: failed to get cooridnates for second point\n";
		return -1;
	    }
	    if (OPS_GetDoubleInput(&ndm, &p3(0)) < 0) {
		opserr<<"WARNING: failed to get cooridnates for third point\n";
		return -1;
	    }
	    if (OPS_GetDoubleInput(&ndm, &p4(0)) < 0) {
		opserr<<"WARNING: failed to get cooridnates for fouth point\n";
		return -1;
	    }

	    // num of particles
	    int numdata = 2;
	    if (OPS_GetIntInput(&numdata, &nump(0)) < 0) {
		opserr<<"WARNING: failed to get particle mesh size\n";
		return -1;
	    }

	} else if (strcmp(geotype,"tri") == 0) {

	    if (OPS_GetNumRemainingInputArgs() < 3*ndm+2) {
		opserr<<"WARNING: insufficient args\n";
		return -1;
	    }

	    // node coord
	    if (OPS_GetDoubleInput(&ndm, &p1(0)) < 0) {
		opserr<<"WARNING: failed to get cooridnates for first point\n";
		return -1;
	    }
	    if (OPS_GetDoubleInput(&ndm, &p2(0)) < 0) {
		opserr<<"WARNING: failed to get cooridnates for second point\n";
		return -1;
	    }
	    if (OPS_GetDoubleInput(&ndm, &p3(0)) < 0) {
		opserr<<"WARNING: failed to get cooridnates for third point\n";
		return -1;
	    }

	    // num of particles
	    int numdata = 2;
	    if (OPS_GetIntInput(&numdata, &nump(0)) < 0) {
		opserr<<"WARNING: failed to get particle mesh size\n";
		return -1;
	    }

	} else if (strcmp(geotype,"line") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 2*ndm+1) {
		opserr<<"WARNING: insufficient args\n";
		return -1;
	    }

	    // node coord
	    if (OPS_GetDoubleInput(&ndm, &p1(0)) < 0) {
		opserr<<"WARNING: failed to get cooridnates for first point\n";
		return -1;
	    }
	    if (OPS_GetDoubleInput(&ndm, &p2(0)) < 0) {
		opserr<<"WARNING: failed to get cooridnates for second point\n";
		return -1;
	    }

	    // num of particles
	    int numdata = 1;
	    if (OPS_GetIntInput(&numdata, &nump(0)) < 0) {
		opserr<<"WARNING: failed to get mesh size\n";
		return -1;
	    }

	} else if (strcmp(geotype,"point") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < ndm) {
		opserr<<"WARNING: insufficient args\n";
		return -1;
	    }
	    
	    // node coord
	    if (OPS_GetDoubleInput(&ndm, &p1(0)) < 0) {
		opserr<<"WARNING: failed to get cooridnates for point\n";
		return -1;
	    }
	    
	} else {
	    opserr<<"WARNING: unknown geometry type\n";
	    return -1;
	}

	// type
	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr<<"WARNING: insufficient args\n";
	    return -1;
	}
	const char* eletype = OPS_GetString();

	// initial velocity
	Vector vel0(ndm);
	if (strcmp(eletype,"-vel") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < ndm) {
		opserr<<"WARNING: insufficient args\n";
		return -1;
	    }
	    
	    // node vel
	    if (OPS_GetDoubleInput(&ndm, &vel0(0)) < 0) {
		opserr<<"WARNING: failed to get initial velocity\n";
		return -1;
	    }

	    // type
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr<<"WARNING: insufficient args\n";
		return -1;
	    }
	    eletype = OPS_GetString();
	}

	// initial pressure
	double p0 = 0;
	if (strcmp(eletype,"-pressure") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr<<"WARNING: insufficient args\n";
		return -1;
	    }
 	    
	    // node pressure
	    int numdata = 1;
	    if (OPS_GetDoubleInput(&numdata, &p0) < 0) {
		opserr<<"WARNING: failed to get initial pressure\n";
		return -1;
	    }

	    // type
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr<<"WARNING: insufficient args\n";
		return -1;
	    }
	    eletype = OPS_GetString();
	}

        // properties
	Vector data;
	if (OPS_GetNumRemainingInputArgs() > 0) {
	    int numdata = OPS_GetNumRemainingInputArgs();
	    data.resize(numdata);
	    if (OPS_GetDoubleInput(&numdata, &data(0)) < 0) {
		opserr<<"WARNING: failed to get fluid properties\n";
		return -1;
	    }
	}

	// get a new group
	ParticleGroup* group = background->newParticleGroup();
	group->setType(eletype);
	group->setProp(data);

	// generate particles
	if (strcmp(geotype,"quad") == 0) {
	  group->qua_d(p1,p2,p3,p4,nump(0),nump(1),vel0,p0);
	} else if (strcmp(geotype,"tri") == 0) {
	    group->tri(p1,p2,p3,nump(0),nump(1),vel0,p0);
	} else if (strcmp(geotype,"line") == 0) {
	    group->line(p1,p2,nump(0),vel0,p0);
	} else if (strcmp(geotype,"point") == 0) {
	    group->point(p1,vel0,p0);
	}

    } else if (strcmp(type,"structure") == 0) {

	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr<<"WARNING: no structural node tags are given\n";
	    return -1;
	}

	int numdata = OPS_GetNumRemainingInputArgs();
	ID tags(numdata);
	if (OPS_GetIntInput(&numdata, &tags(0)) < 0) {
	    opserr<<"WARNING: failed to read structural node tags\n";
	    return -1;
	}
	for (int i=0; i<numdata; i++) {
	    background->addStructuralNode(tags(i));
	}
	
	
    } else if (strcmp(type,"save") == 0) {

	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr<<"WARNING: no file name is given\n";
	    return -1;
	}

	const char* filename = OPS_GetString();
	background->save(filename);

    } else {
	opserr<<"WARNING: unknown background command \n";
	return -1;
    }

    return 0;
}


BackgroundMesh::BackgroundMesh()
    :groups(), fixInfo(), grids(), structuralNodes(), structuralCoord(),
     connectedNodes(), theRecorder(0)
{
}

BackgroundMesh::~BackgroundMesh()
{

    for (int i=0; i<(int)groups.size(); i++) {
    	ParticleGroup* group = groups[i];
    	if (group != 0) delete group;
    }
    groups.clear();

    if (theRecorder != 0) delete theRecorder;
}

ParticleGroup*
BackgroundMesh::newParticleGroup()
{
    ParticleGroup* group = new ParticleGroup;
    groups.push_back(group);
    return group;
}

int
BackgroundMesh::mapParticleToBack()
{
    
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    if (grids.getSize() <= 0) return 0;

    // particles in all grids
    if (particlesInGrids() < 0) {
	opserr<<"WARNING: failed to create grids\n";
	return -1;
    }

    // create nodes in grids which contain particles
    if (nodesInGrids() < 0) {
	opserr<<"WARNING: failed to create background nodes\n";
	return -1;
    }

    // add structural nodes in grids
    if (structureInGrids() < 0) {
    	opserr<<"WARNING: failed to add structure to background mesh\n";
    	return -1;
    }

    // mesh based on grids
    if (mesh() < 0) {
	opserr<<"WARNING: failed to create background elements\n";
	return -1;
    }

    // remove empty elements (no particles in it)
    if (removeEmptyElements() < 0) {
	opserr << "WARNING: failed to remove empty elements\n";
	return -1;
    }

    // fix nodes based on the fixinfo
    if (fix() < 0) {
	opserr<<"WARNING: failed to apply background constraints\n";
	return -1;
    }
    
    return 0;
}

int
BackgroundMesh::mapBackToParticle()
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    if (grids.getSize() <= 0) return 0;
    
    // structure back to grid
    // if (structureToGrids() < 0) {
    // 	opserr<<"WARNING: failed to map structure to background mesh\n";
    // 	return -1;
    // }
    
    // move all the particles
    if (moveParticles() < 0) {
	opserr<<"WARNING: failed to move particles\n";
	return -1;
    }

    // clear all background data
    this->clear();

    return 0;
}

void
BackgroundMesh::addFixInfo(const Vector& min, const Vector& max, const ID& fix)
{
    fixInfo.addInfo(min,max,fix);
}

void
BackgroundMesh::addStructuralNode(int tag)
{
    structuralNodes.insert(tag);
}

int
BackgroundMesh::fix()
{
    
    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
	return 0;
    }

    // for each grid
    grids.reset();
    while (grids.isEnd() == false) {

	// get grid index
	GridIndex index = grids.getIndex();

	// get grid node
	Node* node = grids.getNode(index);

	if (node == 0) {
	    grids.next();
	    continue;
	}
	
	// if its a structure node
	if (structuralNodes.getLocationOrdered(node->getTag()) >= 0) {
	    grids.next();
	    continue;
	}

	// try to fix the grid node
	if (fixInfo.tryFix(node->getTag(),*domain) < 0) {
	    return -1;
	}
	grids.next();
    }

    return 0;
}


int
BackgroundMesh::particlesInGrids()
{
    
    for (int i=0; i<(int)groups.size(); i++) {
	ParticleGroup* group = groups[i];
	if (group == 0) continue;

	for (int j=0; j<group->numParticles(); j++) {
	    Particle* p = group->getParticle(j);
	    if (p == 0) return -1;
	    
	    // locate the particle 
	    const Vector& crds = p->getCrds();
	    if (crds.Size() < 2) return -1;

	    int nx = (int)floor(crds(0)/grids.getSize());
	    int ny = (int)floor(crds(1)/grids.getSize());

	    // add particle to the grid
	    GridIndex index(nx,ny);
	    grids.addParticle(index,p);

	}
    }
    
    return 0;
}



int
BackgroundMesh::structureInGrids()
{
    // get domain
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    
    // mesh size
    double meshsize = this->getMeshsize();
    double halfmeshsize = meshsize/2.0;

    // structural elements
    connectedNodes.clear();
    ElementIter& eles = theDomain->getElements();
    Element* ele = 0;
    while ((ele = eles()) != 0) {
	const ID& elenodes = ele->getExternalNodes();
	for (int i=0; i<elenodes.Size(); i++) {
	    for (int j=0; j<elenodes.Size(); j++) {
		if (j != i) {
		    connectedNodes[elenodes(i)].insert(elenodes(j));
		}
	    }
	}
    }

    // structure nodes
    structuralCoord.clear();
    for (int i=0; i<structuralNodes.Size(); i++) {
	int nodeTag = structuralNodes(i);

	// get node
	Node* snode = theDomain->getNode(nodeTag);
	if (snode == 0) {
	    opserr << "WARNING: node "<<nodeTag<<" does not exist\n";
	    return -1;
	}

	// locate the node 
	const Vector& crds = snode->getCrds();
	const Vector& disp = snode->getTrialDisp();
	if (crds.Size() < 2) {
	    opserr << "WARNING node "<<nodeTag<<" crds size < 2\n";
	    return -1;
	}
	if (disp.Size() < 2) {
	    opserr << "WARNING node "<<nodeTag<<" disp size < 2\n";
	    return -1;
	}

	Vector pos(2);
	pos(0) = crds(0) + disp(0);
	pos(1) = crds(1) + disp(1);

	// store the position
	structuralCoord.push_back(pos);
	
	int nx = (int)floor(pos(0)/grids.getSize());
	int ny = (int)floor(pos(1)/grids.getSize());

	// get GridIndex
	GridIndex index[4];
	index[0] = GridIndex(nx,ny);
	index[1] = index[0].north();
	index[2] = index[0].east();
	index[3] = index[0].northEast();

	// check if any particles in it
	// std::vector<Particle*>* particles = grids.getParticles(index);
	// if (particles == 0 || particles->empty()) {
	//     continue;
	// }

	// the grid node to be replaced
	GridIndex index2;

	double x = index[0].getX(this->getMeshsize());
	double y = index[0].getY(this->getMeshsize());

	// check which grid node the structural node is to replace
	if (pos(0) < x+halfmeshsize) {

	    if (pos(1) < y+halfmeshsize) {
		index2 = index[0];
	    } else {
		index2 = index[1];
	    }
		
	} else {
		
	    if (pos(1) < y+halfmeshsize) {
		index2 = index[2];
	    } else {
		index2 = index[3];
	    }
		
	}

	// check if that grid node exists
	if (grids.hasGrid(index2) == false) {
	    continue;
	}

	// get node of the grid
	Node* gridNode = grids.getNode(index2);
	if (gridNode == 0) {
	    opserr << "WARNING grid exists, but no grid node, that's not normal\n";
	    return -1;
	}

	// check if the grid node is already a structural node
	if (structuralNodes.getLocationOrdered(gridNode->getTag()) >= 0) {

	    // check the structural node is close
	    // to which of the other three grid nodes
	    double dist = -1.0;
	    for (int j=0; j<4; j++) {
		Node* gnode = grids.getNode(index[j]);

		// if the grid node is not a structural node
		if (gnode != 0 && structuralNodes.getLocationOrdered(gnode->getTag())<0) {

		    const Vector& crds1 = gnode->getCrds();
		    const Vector& disp1 = gnode->getDisp();

		    if (crds1.Size() < 2) {
			opserr << "WARNING node "<<gnode->getTag()<<" crds size < 2\n";
			return -1;
		    }

		    if (disp1.Size() < 2) {
			opserr << "WARNING node "<<gnode->getTag()<<" disp size < 2\n";
			return -1;
		    }

		    double dx = pos(0)-crds1(0)-disp1(0);
		    double dy = pos(1)-crds1(1)-disp1(1);
		    double d = dx*dx+dy*dy;
		    if (dist<0 || dist>d) {
			dist = d;
			gridNode = gnode;
			index2 = index[j];
		    }
		}
	    }
	    if (dist < 0) {
		continue;
	    }

	}

	// remove the grid node from domain
	Node* node = theDomain->removeNode(gridNode->getTag());

	// remove the pc of the grid node
	Pressure_Constraint* pc = theDomain->removePressure_Constraint(gridNode->getTag());
	if (pc != 0) delete pc;
	if (node != 0) delete node;

	// replace the grid node with structural node
	grids.setNode(index2, snode);
    }

    return 0;
}


int
BackgroundMesh::nodesInGrids()
{
    
    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
	return 0;
    }

    // for each grid
    grids.reset();
    while (grids.isEnd() == false) {

	// get locations of neibors
	GridIndex index[4];
	index[0] = grids.getIndex();
	if (index[0].isValid() == false) {
	    opserr<<"index[0].isValid() == false\n";
	    grids.next();
	    continue;
	}
	index[1] = index[0].west();
	index[2] = index[0].southWest();
	index[3] = index[0].south();

	// grid's crds
	double x = index[0].getX(grids.getSize());
	double y = index[0].getY(grids.getSize());

	// nodal data
	double wt = 0.0, pressure = 0.0;
	Vector vel;
	
	// map all particles in neighbor to current location
	for (int i=0; i<4; i++) {

	    // particles
	    std::vector<Particle*>* particles = grids.getParticles(index[i]);
	    if (particles == 0) continue;

	    // for each particle
	    for (int j=0; j<(int)particles->size(); j++) {

		Particle* p = (*particles)[j];
		if (p == 0) continue;

		// particle crds
		const Vector& crds = p->getCrds();
		if (crds.Size() < 2) continue;

		// distance from particle to current location
		double dx = crds(0) - x;
		double dy = crds(1) - y;
		double q = sqrt(dx*dx+dy*dy)/grids.getSize();

		// weight for the particle
		double w = QuinticKernel(q, grids.getSize(), crds.Size());

		// add weight
		wt += w;

		// add pressure
		pressure += p->getPressure() * w;

		// add velocity
		const Vector& pvel = p->getVel();
		if (vel.Size() == 0) {
		    vel.resize(pvel.Size());
		    vel.Zero();
		}
		for (int k=0; k<vel.Size(); k++) {
		    if (k < pvel.Size()) {
			vel(k) += w*pvel(k);
		    }
		}
	    }
	}

	if (wt == 0) {
	    grids.next();
	    continue;
	}

	// get nodal data
	pressure /= wt;
	vel /= wt;

	// create a new node
	Node* node = 0;
	node = new Node(findNodeTag(),vel.Size(),x,y);
	
	if (node == 0) {
	    opserr<<"WARNING: run out of memory\n";
	    return -1;
	}
	node->setTrialVel(vel);
	node->commitState();

	// add to domain
	if (domain->addNode(node) == false) {
	    opserr<<"WARNING: failed to add node to domain\n";
	    delete node;
	    return -1;
	}

	// set pressure 
	if (pressure != 0) {
	    Pressure_Constraint* thePC = domain->getPressure_Constraint(node->getTag());
	    if(thePC != 0) {
		thePC->setDomain(domain);
	    } else {
		thePC = new Pressure_Constraint(node->getTag(), 1);
		if(thePC == 0) {
		    opserr<<"WARNING: no enough memory for Pressure_Constraint\n";
		    return -1;
		}
		if(domain->addPressure_Constraint(thePC) == false) {
		    opserr<<"WARNING: failed to add Pressure_Constraint to domain -- ";
		    delete thePC;
		    return -1;
		}
	    }
	    thePC->setPressure(pressure);
	}

	// add node to the grid
	grids.setNode(node);

	// go to next grid
	grids.next();
    }
   
    
    return 0;
}


int
BackgroundMesh::mesh()
{
    
    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
	return 0;
    }

     // for each grid
    grids.reset();
    while (grids.isEnd() == false) {

	// get particles of the grid
	GridIndex index[4];
	index[0] = grids.getIndex();
	if (index[0].isValid() == false) {
	    grids.next();
	    continue;
	}
	std::vector<Particle*>* particles = grids.getParticles();
	if (particles == 0) {
	    grids.next();
	    continue;
	}
	if (particles->empty()) {
	    grids.next();
	    continue;
	}

	// get all four grid index
	index[1] = index[0].east();
	index[2] = index[0].northEast();
	index[3] = index[0].north();

	// get all four nodes
	Node* nodes[4] = {0,0,0,0};
	for (int i=0; i<4; i++) {

	    // get node
	    nodes[i] = grids.getNode(index[i]);
	    if (nodes[i] == 0) return -1;
	}

	// mesh current grid
	ID order;
	meshGrid(nodes, index, order);
	reorder(nodes, order);

	// element type
	ParticleGroup* group = (*particles)[0]->getGroup();
	if (group == 0) {
	    grids.next();
	    continue;
	}
	const char* type = group->getType();
	
	// create element
	if (strcmp(type,"PFEMElement2DBubble") == 0) {
	    const Vector& prop = group->getProp();
	    if (prop.Size() != 6) {
		opserr<<"WARNING: PFEMElement2DBubble element needs 6 args\n";
		return -1;
	    }

	    // add element to domain
	    Element* ele = new PFEMElement2DBubble(findEleTag(),
						   nodes[order(0)]->getTag(),
						   nodes[order(1)]->getTag(),
						   nodes[order(2)]->getTag(),
						   prop(0),prop(1),prop(2),
						   prop(3),prop(4),prop(5));
	    if (ele == 0) {
		opserr<<"WARNING: run out of memory\n";
		return -1;
	    }

	    if (domain->addElement(ele) == false) {
		opserr<<"WARNING: failed to add element to domain\n";
		delete ele;
		return -1;
	    }

	    // add element to background mesh
	    grids.addElement(ele);

	    // check if there is another one
	    if (order.Size() < 6) {
		grids.next();
		continue;
	    }

	    // add element to domain
	    ele = new PFEMElement2DBubble(findEleTag(),
					  nodes[order(3)]->getTag(),
					  nodes[order(4)]->getTag(),
					  nodes[order(5)]->getTag(),
					  prop(0),prop(1),prop(2),
					  prop(3),prop(4),prop(5));
	    if (ele == 0) {
		opserr<<"WARNING: run out of memory\n";
		return -1;
	    }

	    if (domain->addElement(ele) == false) {
		opserr<<"WARNING: failed to add element to domain\n";
		delete ele;
		return -1;
	    }

	    // add element to background mesh
	    grids.addElement(ele);
	}

	// next grid
	grids.next();
    }
    
    return 0;
}

int
BackgroundMesh::meshGrid(Node** nodes, GridIndex* index, ID& order)
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    
    order.resize(6);

    // element 012 and 023
    order(0) = 0;
    order(1) = 1;
    order(2) = 2;
    order(3) = 0;
    order(4) = 2;
    order(5) = 3;

    const Vector& crds0 = nodes[0]->getCrds();
    const Vector& crds1 = nodes[1]->getCrds();
    const Vector& crds2 = nodes[2]->getCrds();
    const Vector& crds3 = nodes[3]->getCrds();

    // check how many structural node
    bool structure[4] = {false, false, false, false};
    int nums = 0;
    for (int i=0; i<4; i++) {
	if (structuralNodes.getLocationOrdered(nodes[i]->getTag()) >= 0) {
	    structure[i] = true;
	    nums += 1;
	}
    }

    // if there are 3 structural nodes
    if (nums == 3) {

	if (structure[0] == false) {
	    
	    // if 1 and 3 are connected
	    if (connectedNodes[nodes[1]->getTag()].getLocationOrdered(nodes[3]->getTag()) >= 0) {
		order.resize(3);
		order(0) = 0;
		order(1) = 1;
		order(2) = 3;
	    } else {
		// keep same
	    }
	    return 0;
	    
	} else if (structure[1] == false) {

	    // if 0 and 2 are connected
	    if (connectedNodes[nodes[0]->getTag()].getLocationOrdered(nodes[2]->getTag()) >= 0) {
		order.resize(3);
		order(0) = 0;
		order(1) = 1;
		order(2) = 2;
	    } else {
		order(2) = 3;
		order(3) = 1;
	    }
	    return 0;

	} else if (structure[2] == false) {

	    // if 1 and 3 are connected
	    if (connectedNodes[nodes[1]->getTag()].getLocationOrdered(nodes[3]->getTag()) >= 0) {
		order.resize(3);
		order(0) = 1;
		order(1) = 2;
		order(2) = 3;
	    } else {
		// keep same
	    }
	    return 0;
	    
	} else if (structure[3] == false) {

	    // if 0 and 2 are connected
	    if (connectedNodes[nodes[0]->getTag()].getLocationOrdered(nodes[2]->getTag()) >= 0) {
		order.resize(3);
		order(0) = 0;
		order(1) = 2;
		order(2) = 3;
	    } else {
		order(2) = 3;
		order(3) = 1;
	    }
	    return 0;
	}
	
    }

    // if 1 and 3 are structural nodes
    if (nums == 2 && structure[1] && structure[3]) {

	// if 1 and 3 connected
	if (connectedNodes[nodes[1]->getTag()].getLocationOrdered(nodes[3]->getTag()) >= 0) {
	    order(2) = 3;
	    order(3) = 1;
	    return 0;
	}
    } else if (nums == 2 && structure[0] && structure[2]) {
	
	// if 0 and 2 connected
	if (connectedNodes[nodes[0]->getTag()].getLocationOrdered(nodes[2]->getTag()) >= 0) {
	    return 0;
	}
    }

    // for all else having structure, select the shorter diagonal
    if (nums > 0) {
	Vector v02 = crds0-crds2;
	Vector v13 = crds1-crds3;

	double d02 = 0.0;
	for (int i=0; i<v02.Size(); i++) {
	    d02 += v02(i)*v02(i);
	}

	double d13 = 0.0;
	for (int i=0; i<v13.Size(); i++) {
	    d13 += v13(i)*v13(i);
	}

	if (d13 < d02) {
	    order(2) = 3;
	    order(3) = 1;
	}

	return 0;
	
    }

    // for all fluid,
    // if node 1 or 3 is corner, change to element 013 and 123
    if (grids.isCorner(index[1]) || grids.isCorner(index[3])) {
	order(2) = 3;
	order(3) = 1;
    }

    return 0;
}

int
BackgroundMesh::reorder(Node** nodes, ID& order)
{
    if (order.Size() < 3) return 0;

    double x[3], y[3];
    for (int i=0; i<3; i++) {
	const Vector& crds = nodes[order(i)]->getCrds();
	const Vector& disp = nodes[order(i)]->getDisp();
	x[i] = crds(0) + disp(0);
	y[i] = crds(1) + disp(1);
    }

    if (x[1]*y[2]-x[2]*y[1]+x[2]*y[0]-x[0]*y[2]+x[0]*y[1]-x[1]*y[0] < 0) {
	int temp = order(1);
	order(1) = order(2);
	order(2) = temp;
    }

    if (order.Size() < 6) return 0;

    for (int i=0; i<3; i++) {
	const Vector& crds = nodes[order(i+3)]->getCrds();
	const Vector& disp = nodes[order(i+3)]->getDisp();
	x[i] = crds(0) + disp(0);
	y[i] = crds(1) + disp(1);
    }

    if (x[1]*y[2]-x[2]*y[1]+x[2]*y[0]-x[0]*y[2]+x[0]*y[1]-x[1]*y[0] < 0) {
	int temp = order(4);
	order(4) = order(5);
	order(5) = temp;
    }

    return 0;
}

int
BackgroundMesh::removeEmptyElements()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    
    // for each grid
    grids.reset();
    while (grids.isEnd() == false) {

	// get grid elements
	std::vector<Element*>* eles = grids.getElements();
	if (eles == 0 || eles->empty()) {
	    grids.next();
	    continue;
	}

	// get grid four corners
	GridIndex index[4];
	index[0] = grids.getIndex();
	if (index[0].isValid() == false) {
	    grids.next();
	    continue;
	}
	index[1] = index[0].east();
	index[2] = index[0].northEast();
	index[3] = index[0].north();

	// get grid four nodes
	Node* nodes[4];
	for (int i=0; i<4; i++) {
	    nodes[i] = grids.getNode(index[i]);
	    if (nodes[i] == 0) {
		opserr << "WARNING: no corner node \n";
		return -1;
	    }
	}

	// if a corner node is structural node
	// gather particles
	std::map< GridIndex, std::vector<Particle*>* > particles;
	for (int i=0; i<4; i++) {
	    if (structuralNodes.getLocationOrdered(nodes[i]->getTag()) >= 0) {
		particles[index[i]] = grids.getParticles(index[i]);
		particles[index[i].west()] = grids.getParticles(index[i].west());
		particles[index[i].southWest()] = grids.getParticles(index[i].southWest());
		particles[index[i].south()] = grids.getParticles(index[i].south());
	    }
	}
	if (particles.empty()) {
	    grids.next();
	    continue;
	}

	// check if any grid elements are empty
	Vector N;
	for (int i=0; i<(int)eles->size(); i++) {

	    // get element
	    Element* ele = (*eles)[i];
	    if (ele == 0) {
		continue;
	    }
	    if (ele->getNumExternalNodes() != 6) {
		opserr << "WARNING: fluid element has number of fluid nodes != 3\n";
		return -1;
	    }

	    // get element nodes 
	    Node** elenodes = ele->getNodePtrs();
	    const Vector& crds1 = elenodes[0]->getCrds();
	    const Vector& crds2 = elenodes[2]->getCrds();
	    const Vector& crds3 = elenodes[4]->getCrds();
	    const Vector& disp1 = elenodes[0]->getDisp();
	    const Vector& disp2 = elenodes[2]->getDisp();
	    const Vector& disp3 = elenodes[4]->getDisp();

	    if (crds1.Size()<2 || crds2.Size()<2 || crds3.Size()<2 ||
		disp1.Size()<2 || disp2.Size()<2 || disp3.Size()<2) {
		opserr << "WARNING: dim < 2 or ndf < 2\n";
		return -1;
	    }

	    double x1 = crds1(0)+disp1(0);
	    double y1 = crds1(1)+disp1(1);
	    double x2 = crds2(0)+disp2(0);
	    double y2 = crds2(1)+disp2(1);
	    double x3 = crds3(0)+disp3(0);
	    double y3 = crds3(1)+disp3(1);

	    // check each particle
	    bool empty = true;
	    std::map< GridIndex, std::vector<Particle*>* >::iterator it;
	    for (it=particles.begin(); it!=particles.end(); it++) {
		if (it->second == 0) continue;
		if (it->second->empty()) continue;

		for (int j=0; j<(int)it->second->size(); j++) {
		    Particle* p = (*(it->second))[j];
		    if (p == 0) continue;
		    const Vector& crds = p->getCrds();
		    if (crds.Size() < 2) continue;
		    getNForTri(x1,y1,x2,y2,x3,y3,crds(0),crds(1),N);

		    // this particle is in the element
		    if (N(0)>=0 && N(0)<=1 && N(1)>=0 && N(1)<=1 && N(2)>=0 && N(2)<=1) {
			empty = false;
			break;
		    }
		}
	    }

	    // if empty element, remove it
	    if (empty) {
		theDomain->removeElement(ele->getTag());
		delete ele;
		(*eles)[i] = 0;
	    }
	    
	}
	
	grids.next();
    }
    return 0;
}

double
BackgroundMesh::QuinticKernel(double q, double h, int ndm)
{
    if (q<0 || q>2) return 0.0;
    double aD = 0.0;
    if (ndm == 2) {
	aD = 7.0/(4*pi*h*h);
    } else if (ndm == 3) {
	aD = 7.0/(8*pi*h*h*h);
    }
    double a = 1.0-q/2.0;

    return aD*a*a*a*a*(2*q+1);
}

int
BackgroundMesh::findNodeTag()
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    
    // get node tag
    NodeIter& theNodes = domain->getNodes();
    Node* theNode = theNodes();
    int currtag = 0;
    if(theNode != 0) currtag = theNode->getTag();

    return currtag-1;
}

int
BackgroundMesh::findEleTag()
{
    
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    
    // get ele tag
    ElementIter& theEles = domain->getElements();
    Element* theEle = theEles();
    int currtag = 0;
    if(theEle != 0) currtag = theEle->getTag();

    return currtag-1;
}

// int
// BackgroundMesh::structureToGrids()
// {
    
    // Domain* theDomain = OPS_GetDomain();
    // if (theDomain == 0) return -1;
    
    // structure.reset();
    // while(structure.isEnd() == false) {

    // 	// grid and structure node
    // 	GridIndex index = structure.getGridIndex();
    // 	Node* nodes[9] = {0,0,0,0,0,0,0,0,0};
    // 	nodes[8] = structure.getNode();

    // 	if (nodes[8] == 0) {
    // 	    opserr<<"WARNING: structural node not exists\n";
    // 	    return -1;
    // 	}

    // 	// neighbor grid nodes
    // 	nodes[0] = grids.getNode(index.northEast());
    // 	nodes[1] = grids.getNode(index.north());
    // 	nodes[2] = grids.getNode(index.northWest());
    // 	nodes[3] = grids.getNode(index.west());
    // 	nodes[4] = grids.getNode(index.southWest());
    // 	nodes[5] = grids.getNode(index.south());
    // 	nodes[6] = grids.getNode(index.southEast());
    // 	nodes[7] = grids.getNode(index.east());

    // 	// grid node data
    // 	double wt = 0.0, gp = 0.0;
    // 	Vector gvel;

    // 	// grid coordinates
    // 	double x = index.getX(this->getMeshsize());
    // 	double y = index.getY(this->getMeshsize());

    // 	// number of neighbor nodes
    // 	int numNeigh = 0;
    // 	for (int i=0; i < 9; i++) {
    // 	    if (nodes[i] != 0) {
    // 		numNeigh++;
    // 	    }
    // 	}

    // 	// map neighbor grid nodes and structural node to the grid
    // 	for (int i=0; i < 9; i++) {
    // 	    if (nodes[i] == 0) continue;

    // 	    // neighbor node coordinates
    // 	    const Vector& crds = nodes[i]->getCrds();
    // 	    if (crds.Size() < 2) continue;

    // 	    double nex = crds(0);
    // 	    double ney = crds(1);

    // 	    // structural node
    // 	    if (i == 8) {
    // 		const Vector& disp = nodes[i]->getTrialDisp();
    // 		if (disp.Size() < 2) {
    // 		    continue;
    // 		}
    // 		nex += disp(0);
    // 		ney += disp(1);
    // 	    }

    // 	    // distance from neighbor node to current grid node
    // 	    double dx = nex - x;
    // 	    double dy = ney - y;
    // 	    double q = sqrt(dx*dx+dy*dy)/numNeigh;

    // 	    // weight for the the neighbor node
    // 	    double w = QuinticKernel(q, numNeigh, crds.Size());

    // 	    // add weight
    // 	    wt += w;

    // 	    // add pressure
    // 	    Pressure_Constraint* pc = theDomain->getPressure_Constraint(nodes[i]->getTag());
    // 	    if (pc != 0) {
    // 		gp += pc->getPressure() * w;
    // 	    }

    // 	    // add velocity
    // 	    const Vector& vel = nodes[i]->getTrialVel();
    // 	    if (gvel.Size() == 0) {
    // 		gvel.resize(vel.Size());
    // 		gvel.Zero();
    // 	    }
    // 	    for (int k=0; k<gvel.Size(); k++) {
    // 		if (k < vel.Size()) {
    // 		    gvel(k) += w*vel(k);
    // 		}
    // 	    }
    // 	}

    // 	// get nodal data
    // 	if (wt == 0) {
    // 	    gp = 0.0;
    // 	    gvel.Zero();
    // 	}
    // 	gp /= wt;
    // 	gvel /= wt;
	    
    // 	// set the grid node
    // 	Node* gnode = new Node(findNodeTag(), gvel.Size(), x, y);
    // 	if (gnode == 0) {
    // 	    opserr<<"WARNING: failed to create grid node run out of memory\n";
    // 	    return -1;
    // 	}
    // 	gnode->setTrialVel(gvel);
    // 	gnode->commitState();

    // 	if (theDomain->addNode(gnode) == false) {
    // 	    opserr << "WARNING: failed to add node to domain\n";
    // 	    delete gnode;
    // 	    return -1;
    // 	}
	
    // 	// set the pc for grid node
    // 	Pressure_Constraint* pc = theDomain->getPressure_Constraint(gnode->getTag());
    // 	// opserr << "grid node "<<gnode->getTag()<<" ndf = "<<gnode->getNumberDOF()<<"\n";
    // 	// Node* theNode = theDomain->getNode(gnode->getTag());
    // 	// if (theNode != 0) {
    // 	//     opserr << "grid node "<<gnode->getTag()<<" exists\n";
    // 	// } else {
    // 	//     opserr << "grid node "<<gnode->getTag()<<" not exists\n";
    // 	// }
    // 	if (pc == 0) {
    // 	    pc = new Pressure_Constraint(gnode->getTag(), 1);
    // 	    if(pc == 0) {
    // 		opserr<<"WARNING: no enough memory for Pressure_Constraint\n";
    // 		return -1;
    // 	    }
    // 	    if(theDomain->addPressure_Constraint(pc) == false) {
    // 		opserr<<"WARNING: failed to add Pressure_Constraint to domain -- ";
    // 		delete pc;
    // 		return -1;
    // 	    }
    // 	}
    // 	pc->setPressure(gp);
	
    // 	// switch back grid and structural nodes
    // 	grids.setNode(index, gnode);

    // 	structure.next();
    // }
    
//     return 0;
// }

int
BackgroundMesh::moveParticles()
{
    
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    
    // for each grid
    grids.reset();
    while (grids.isEnd() == false) {

	// get all four grid points
	GridIndex index[4];
	index[0] = grids.getIndex();
	if (index[0].isValid() == false) {
	    grids.next();
	    continue;
	}
	index[1] = index[0].east();
	index[2] = index[0].northEast();
	index[3] = index[0].north();

		    
	// get grid elements
	std::vector<Element*>* eles = grids.getElements();
	if (eles == 0 || eles->empty()) {
	    grids.next();
	    continue;
	}


	// get all four nodes
	Node* nodes[4];
	for (int i=0; i<4; i++) {

	    // get node
	    nodes[i] = grids.getNode(index[i]);
	    if (nodes[i] == 0) {
		opserr << "WARNING: no corner node "<<i<<"\n";
		return -1;
	    }
	}

	// if a corner node is structural node
	// gather particles
	std::map< GridIndex, std::vector<Particle*>* > allparticles;
	for (int i=0; i<4; i++) {
	    if (structuralNodes.getLocationOrdered(nodes[i]->getTag()) >= 0) {
		allparticles[index[i]] = grids.getParticles(index[i]);
		allparticles[index[i].west()] = grids.getParticles(index[i].west());
		allparticles[index[i].southWest()] = grids.getParticles(index[i].southWest());
		allparticles[index[i].south()] = grids.getParticles(index[i].south());
	    }
	}

	// if no structural node, move particles in the cell
	if (allparticles.empty()) {

	    std::vector<Particle*>* particles = grids.getParticles();
	    if (particles == 0) {
		grids.next();
		continue;
	    }
	    if (particles->empty()) {
		grids.next();
		continue;
	    }

	    // grid's crds
	    double x0 = index[0].getX(grids.getSize());
	    double y0 = index[0].getY(grids.getSize());

	    // shape functions
	    Vector N;
	
	    // move all particles in the cell
	    for (int i=0; i<(int)particles->size(); i++) {

		// get shape function
		const Vector& crds = (*particles)[i]->getCrds();
		if (crds.Size() < 2) continue;
		getNForRect(x0,y0,grids.getSize(),grids.getSize(),crds(0),crds(1),N);
		Vector pdisp(crds.Size());
		Vector pvel((*particles)[i]->getVel().Size());
		double ppre = 0.0;

		// interpolation
		for (int j=0; j<4; j++) {
		    const Vector& disp = nodes[j]->getDisp();
		    const Vector& vel = nodes[j]->getVel();
		    Pressure_Constraint* pc = domain->getPressure_Constraint(nodes[j]->getTag());
		    double pressure = 0.0;
		    if (pc != 0) pressure = pc->getPressure();

		    for (int k=0; k<pdisp.Size(); k++) {
			if (k < disp.Size()) {
			    pdisp(k) += N(j)*disp(k);
			}
		    }
		    for (int k=0; k<pvel.Size(); k++) {
			if (k < vel.Size()) {
			    pvel(k) += N(j)*vel(k);
			}
		    }
		    if (pressure != 0.0) {
			ppre += N(j) * pressure;
		    }
		}

		// move the particle
		(*particles)[i]->move(pdisp);
		(*particles)[i]->setVel(pvel);
		(*particles)[i]->setPressure(ppre);
	    }
	} else {

	    // if there is structural nodes, move particles in elements

	    // move particles in the elements
	    Vector N;
	    for (int i=0; i<(int)eles->size(); i++) {

		// get element
		Element* ele = (*eles)[i];
		if (ele == 0) continue;
		if (ele->getNumExternalNodes() != 6) {
		    opserr << "WARNING: fluid element has number of fluid nodes != 3\n";
		    return -1;
		}

		// get element nodes original coordinates
		Node** elenodes = ele->getNodePtrs();
		double pressure[3] = {0.,0.,0.};
		double x[3], y[3], vx[3], vy[3];
		double posx[3], posy[3];
		for (int j=0; j<3; j++) {
		    int loc = structuralNodes.getLocationOrdered(elenodes[2*j]->getTag());

                    // get displacement
		    const Vector& disp = elenodes[2*j]->getDisp();
		    double ux = disp(0);
		    double uy = disp(1);

		    // node cooridnates of original step
		    const Vector& crds = elenodes[2*j]->getCrds();
		    x[j] = crds(0);
		    y[j] = crds(1);

		    // current position for both fluid and structure
		    posx[j] = x[j] + ux;
		    posy[j] = y[j] + uy;

		    // structural node coordinates of last time step
		    if (loc >= 0) {
			const Vector& crds1 = structuralCoord[loc];
			x[j] = crds1(0);
			y[j] = crds1(1);
		    }

		    // get vel
		    const Vector& vel = elenodes[2*j]->getVel();
		    vx[j] = vel(0);
		    vy[j] = vel(1);

		    // get pressure
		    Pressure_Constraint* pc = domain->getPressure_Constraint(elenodes[2*j]->getTag());
		    if (pc != 0) {
			pressure[j] = pc->getPressure();
		    }
		    
		}

		// check each particle
		std::map< GridIndex, std::vector<Particle*>* >::iterator it;
		for (it=allparticles.begin(); it!=allparticles.end(); it++) {
		    if (it->second == 0) continue;
		    if (it->second->empty()) continue;

		    for (int j=0; j<(int)it->second->size(); j++) {
			Particle* p = (*(it->second))[j];
			if (p == 0) continue;
			const Vector& crds = p->getCrds();
			if (crds.Size() < 2) continue;
			getNForTri(x[0],y[0],x[1],y[1],x[2],y[2],crds(0),crds(1),N);

			// this particle is in the element
			if (N(0)>=0 && N(0)<=1 && N(1)>=0 && N(1)<=1 && N(2)>=0 && N(2)<=1) {

			    // move the particle
			    Vector pcrds(crds.Size());
			    Vector pvel(p->getVel().Size());
			    double ppre = 0.0;

			    // interpolation
			    for (int k=0; k<3; k++) {

				pcrds(0) += N(k) * posx[k];
				pcrds(1) += N(k) * posy[k];
				pvel(0) += N(k) * vx[k];
				pvel(1) += N(k) * vy[k];
				ppre += N(k) * pressure[k];
			    }

			    // move the particle
			    p->moveTo(pcrds);
			    p->setVel(pvel);
			    p->setPressure(ppre);
			}
		    }
		}
		
	    }

	}
	
	grids.next();
    }
    
    return 0;
}

void
BackgroundMesh::getNForRect(double x0, double y0, double hx, double hy, double x, double y,
			    Vector& N)
{
    // compute local coordinate of the particle
    double xl = (x-x0)/hx;
    double yl = (y-y0)/hy;

    // map to [-1, 1]
    xl = xl*2-1;
    yl = yl*2-1;

    // shape function
    N.resize(4);
    N(0) = (1-xl)*(1-yl)/4.0;
    N(1) = (1+xl)*(1-yl)/4.0;
    N(2) = (1+xl)*(1+yl)/4.0;
    N(3) = (1-xl)*(1+yl)/4.0;
}

void
BackgroundMesh::getNForTri(double x1, double y1, double x2, double y2,
			   double x3, double y3, double x, double y,
			   Vector& N)
{
    double a1 = x2*y3-x3*y2;
    double a2 = x3*y1-x1*y3;
    double a3 = x1*y2-x2*y1;

    double b1 = y2-y3;
    double b2 = y3-y1;
    double b3 = y1-y2;

    double c1 = x3-x2;
    double c2 = x1-x3;
    double c3 = x2-x1;

    double A =a1+a2+a3;

    N.resize(3);
    N(0) = (a1+b1*x+c1*y)/A;
    N(1) = (a2+b2*x+c2*y)/A;
    N(2) = (a3+b3*x+c3*y)/A;
}

void
BackgroundMesh::clear()
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return;

    
    fixInfo.clear(*domain);
    grids.clear(structuralNodes);

    // remove structural pressure nodes 
    for (int i=0; i<structuralNodes.Size(); i++) {
	// remove pc
	Pressure_Constraint* pc = domain->removePressure_Constraint(structuralNodes(i));
	if (pc != 0) delete pc;
    }
}

void
BackgroundMesh::save(const char* filename)
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return;

    // create particle file
    std::string particlefilename(filename);
    particlefilename += ".particle";
    std::ofstream outfile(particlefilename.c_str(), std::fstream::app);

    // get total num of particles
    int numParticles = 0;

    for (int i=0; i<(int)groups.size(); i++) {
	ParticleGroup* group = groups[i];
	if (group == 0) continue;
	int num = group->numParticles();
	numParticles += num;
    }

    // first line: time, meshsize, numParticles
    outfile << domain->getCurrentTime() << " " << grids.getSize() << " ";
    outfile << numParticles << " ";
    outfile << 0 << " " << 0 <<"\n";

    // each group 
    for (int i=0; i<(int)groups.size(); i++) {
	ParticleGroup* group = groups[i];
	if (group == 0) continue;
	
	int num = group->numParticles();

	// each particle
	for (int j=0; j<num; j++) {
	    Particle* p = group->getParticle(j);
	    if (p == 0) {
		opserr << "WARNING: particle dose not exist\n";
		return;
	    }

	    // particle line : x y (z) vx vy (vz) pressure
	    const Vector& crds = p->getCrds();
	    const Vector& vel = p->getVel();
	    double pressure = p->getPressure();

	    for (int k = 0; k < crds.Size(); k++) {
		outfile << crds(k) << " ";
	    }
	    for (int k = 0; k < vel.Size(); k++) {
		outfile << vel(k) << " ";
	    }
	    outfile << pressure << "\n";
	}
    }
    outfile.close();

    // number of elements and nodes
    int numEles = 0;
    std::map<int,int> nodetags;
    grids.reset();
    while (grids.isEnd() == false) {

	// get elements
	std::vector<Element*>* eles = grids.getElements();
	if (eles == 0) {
	    grids.next();
	    continue;
	}

	// add num
	int num = (int)eles->size();
	numEles += num;

	// node tags
	for (int i=0; i<num; i++) {

	    // elemental nodes
	    const ID& elenodes = (*eles)[i]->getExternalNodes();

	    // each ele node
	    for (int j=0; j<elenodes.Size(); j+=2) {

		// insert node tag
		nodetags[elenodes(j)] = 0;
	    }
	}
	grids.next();
    }

    // create node file
    std::string nodefilename(filename);
    nodefilename += ".node";
    outfile.open(nodefilename.c_str(), std::fstream::app);

    // first line: time, numNodes
    outfile << domain->getCurrentTime() << " " << (int)nodetags.size() << " " ;
    outfile << 0 << 0 << 0<<"\n";

    // saving
    int index = 0;
    for (std::map<int,int>::iterator it=nodetags.begin(); it!=nodetags.end(); it++) {

	// set index
	it->second = index++;

	// get node
	Node* node = domain->getNode(it->first);
	if (node == 0) {
	    opserr<<"WARNING: node "<<it->first<<" does not exist\n";
	    return;
	}

	// save nodal data
	const Vector& crds = node->getCrds();
	const Vector& vel = node->getVel();
	if (crds.Size() < 2 || vel.Size() < 2) {
	    opserr<<"WARNING: crds and vel size < 2\n";
	    return;
	}
	double pressure = 0.0;
	Pressure_Constraint* pc = domain->getPressure_Constraint(it->first);
	if (pc != 0) {
	    pressure = pc->getPressure();
	}
	outfile << crds(0) << " " << crds(1) << " ";
	outfile << vel(0) << " " << vel(1) << " ";
	outfile << pressure << "\n";
    }
    outfile.close();
    
    // create element file
    std::string elefilename(filename);
    elefilename += ".ele";
    outfile.open(elefilename.c_str(), std::fstream::app);

    // first line: time, numEles
    outfile << domain->getCurrentTime() << " " << numEles << " " ;
    outfile << 0 <<"\n";

    // saving
    grids.reset();
    while (grids.isEnd() == false) {

	// get elements
	std::vector<Element*>* eles = grids.getElements();
	if (eles == 0) {
	    grids.next();
	    continue;
	}

	// save
	int num = (int)eles->size();
	for (int i=0; i<num; i++) {

	    // elemental nodes
	    const ID& elenodes = (*eles)[i]->getExternalNodes();

	    // each ele node
	    for (int j=0; j<elenodes.Size(); j+=2) {

		// get index
		int index = nodetags[elenodes(j)];
		outfile << index << " ";
	    }
	    outfile << "\n";
	}

	grids.next();
    }
    outfile.close();

}

void
BackgroundMesh::setRecorder(PVDRecorder* recorder)
{
    theRecorder = recorder;
    if (theRecorder != 0) {
	Domain* theDomain = OPS_GetDomain();
	if (theDomain == 0) return;
	theRecorder->setDomain(*theDomain);
    }
}

int
BackgroundMesh::record(int ctag, double timestamp)
{
    if (theRecorder == 0) {
	opserr << "WARNING no pvd recorder is set for background mesh\n";
	return -1;
    }

    return theRecorder->record(ctag, timestamp);
}
