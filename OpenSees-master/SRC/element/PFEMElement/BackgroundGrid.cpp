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

#include "BackgroundGrid.h"
#include <Domain.h>
#include <Element.h>
#include <Node.h>
#include <Pressure_Constraint.h>
#include <elementAPI.h>

BackgroundGrid::BackgroundGrid()
    :data(), size(0), iter()
{
    iter = data.end();
}

BackgroundGrid::~BackgroundGrid()
{
    // for each grid
    reset();
    while (isEnd() == false) {

	// remove data
	if (iter->second != 0) {
	    delete iter->second;
	    iter->second = 0;
	}

	next();
    }
    
    data.clear();
    iter = data.end();
}

void
BackgroundGrid::reset(const GridIndex& index)
{
    iter = data.find(index);
}

void
BackgroundGrid::addParticle(const GridIndex& index, Particle* p)
{
    // create grid
    GridData* griddata = new GridData;

    // insert the grid
    std::pair<std::map<GridIndex,GridData*>::iterator, bool> res;
    res = data.insert(std::make_pair(index, griddata));

    // if the grid already existed
    if (res.second == false) {
	delete griddata;
	griddata = res.first->second;
    }

    // add particle
    griddata->particles.push_back(p);

    // other grids
    griddata = new GridData;
    if (data.insert(std::make_pair(index.east(), griddata)).second == false) {
	delete griddata;
    }

    griddata = new GridData;
    if (data.insert(std::make_pair(index.north(), griddata)).second == false) {
	delete griddata;
    }

    griddata = new GridData;
    if (data.insert(std::make_pair(index.northEast(), griddata)).second == false) {
	delete griddata;
    }
}

void
BackgroundGrid::addElement(const GridIndex& index, Element* e)
{
    data[index]->elements.push_back(e);
}

void
BackgroundGrid::setNode(const GridIndex& index, Node* nd)
{
    data[index]->node = nd;
}

std::vector<Particle*>*
BackgroundGrid::getParticles(const GridIndex& index)
{
    std::map<GridIndex,GridData*>::iterator it = data.find(index);
    if (it == data.end()) return 0;
    
    return &(it->second->particles);
}

std::vector<Element*>*
BackgroundGrid::getElements(const GridIndex& index)
{
    std::map<GridIndex,GridData*>::iterator it = data.find(index);
    if (it == data.end()) return 0;

    return &(it->second->elements);
}

Node*
BackgroundGrid::getNode(const GridIndex& index)
{
    std::map<GridIndex,GridData*>::iterator it = data.find(index);
    if (it == data.end()) return 0;
    return it->second->node;
}

bool
BackgroundGrid::hasGrid(const GridIndex& index)
{
    std::map<GridIndex,GridData*>::iterator it = data.find(index);
    if (it == data.end()) return false;
    return true;
}

void
BackgroundGrid::clear(const ID& structuralNodes)
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return;
    
    // for each grid
    reset();
    while (isEnd() == false) {

	// remove elements
	std::vector<Element*>* eles = getElements();
	    
	if (eles != 0) {
	    for (int i=0; i<(int)eles->size(); i++) {
		Element* ele = (*eles)[i];
		if (ele != 0) {
		    if (theDomain->removeElement(ele->getTag()) != 0) {
			delete ele;
		    }
		}
	    }
	}
	    
	// remove nodes
	Node* node = getNode();
	if (node != 0) {

	    if (structuralNodes.getLocationOrdered(node->getTag()) < 0) {
		// remove pc
		Pressure_Constraint* pc = theDomain->removePressure_Constraint(node->getTag());
		if (pc != 0) delete pc;
		
		if (theDomain->removeNode(node->getTag()) != 0) {
		    delete node;
		}
	    }
	}

	// remove data
	if (iter->second != 0) {
	    delete iter->second;
	    iter->second = 0;
	}

	next();
    }
    
    data.clear();
    iter = data.end();
}

GridIndex
BackgroundGrid::getIndex() const
{
    if (iter == data.end()) return GridIndex();
    return iter->first;
}

std::vector<Particle*>*
BackgroundGrid::getParticles()
{
    if (iter == data.end()) return 0;
    if (iter->second == 0) return 0;
    return &(iter->second->particles);
}

std::vector<Element*>*
BackgroundGrid::getElements()
{
    if (iter == data.end()) return 0;
    if (iter->second == 0) return 0;
    return &(iter->second->elements);
}

Node*
BackgroundGrid::getNode()
{
    if (iter == data.end()) return 0;
    if (iter->second == 0) return 0;
    return iter->second->node;
}

void
BackgroundGrid::addParticle(Particle* p)
{
    if (iter == data.end()) return;
    if (iter->second == 0) return;
    iter->second->particles.push_back(p);
}

void
BackgroundGrid::addElement(Element* e)
{
    if (iter == data.end()) return;
    if (iter->second == 0) return;
    iter->second->elements.push_back(e);
}

void
BackgroundGrid::setNode(Node* nd)
{
    if (iter == data.end()) return;
    if (iter->second == 0) return;
    iter->second->node = nd;
}

void
BackgroundGrid::next()
{
    if (isEnd()) return;
    iter++;
}

bool
BackgroundGrid::isCorner(const GridIndex& center) const
{
    std::map<GridIndex,GridData*>::const_iterator iternorth = data.find(center.north());
    std::map<GridIndex,GridData*>::const_iterator itersouth = data.find(center.south());
    std::map<GridIndex,GridData*>::const_iterator itereast = data.find(center.east());
    std::map<GridIndex,GridData*>::const_iterator iterwest = data.find(center.west());

    if (iternorth == data.end() && iterwest == data.end()) return true;
    if (iternorth == data.end() && itereast == data.end()) return true;
    if (itersouth == data.end() && iterwest == data.end()) return true;
    if (itersouth == data.end() && itereast == data.end()) return true;

    return false;
}
