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

#ifndef BackgroundStructure_h
#define BackgroundStructure_h

#include "BackgroundGrid.h"
#include <ID.h>
#include <Node.h>
#include <map>

class BackgroundStructure
{

public:
    BackgroundStructure();
    ~BackgroundStructure();

    void addNode(int nodeTag);
    void setIndex(const GridIndex& index, int nodeTag);
    const ID& getNodes() const {return nodes;}
    bool hasGrid(const GridIndex& index);

    Node* getNode();
    GridIndex getGridIndex() const;

    void reset();
    void reset(const GridIndex& index);
    bool isEnd() const;
    void next();

    void clear();
    
private:
    ID nodes;
    std::map<GridIndex,int> nodeGridIndex;
    std::map<GridIndex,int>::iterator iter;
};


#endif
