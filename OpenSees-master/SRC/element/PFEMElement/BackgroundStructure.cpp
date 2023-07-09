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

#include "BackgroundStructure.h"
#include <elementAPI.h>
#include <Domain.h>

BackgroundStructure::BackgroundStructure()
    :nodes(), nodeGridIndex(), iter()
{
    iter = nodeGridIndex.end();
}

BackgroundStructure::~BackgroundStructure()
{
    this->clear();
}

void
BackgroundStructure::addNode(int nodeTag)
{
    nodes.insert(nodeTag);
}

void
BackgroundStructure::reset()
{
    iter = nodeGridIndex.begin();
}

void
BackgroundStructure::reset(const GridIndex& index)
{
    iter = nodeGridIndex.find(index);
}

bool
BackgroundStructure::isEnd() const
{
    return iter == nodeGridIndex.end();
}

void
BackgroundStructure::next()
{
    if (isEnd()) return;
    iter++;
}

Node*
BackgroundStructure::getNode()
{
    if (isEnd()) return 0;

    if (iter->second == 0) return 0;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return 0;

    return theDomain->getNode(iter->second);
}

GridIndex
BackgroundStructure::getGridIndex() const
{
    if (isEnd()) return GridIndex();

    if (iter->second == 0) return GridIndex();

    return iter->first;
}

void
BackgroundStructure::setIndex(const GridIndex& index, int nodeTag)
{
    // no two structural nodes replace same grid node
    nodeGridIndex[index] = nodeTag;
}

bool
BackgroundStructure::hasGrid(const GridIndex& index)
{
    std::map<GridIndex,int>::iterator it = nodeGridIndex.find(index);
    if (it == nodeGridIndex.end()) return false;
    return true;
}

void
BackgroundStructure::clear()
{
    nodeGridIndex.clear();
    iter = nodeGridIndex.end();
}
