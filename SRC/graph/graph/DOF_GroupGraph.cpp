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
                                                                        
// $Revision: 1.3 $
// $Date: 2005-11-03 23:11:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/graph/DOF_GroupGraph.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the class definition for DOF_GroupGraph.
// DOF_GroupGraph is a graph of the DOF_Groups in the domain. It is used by the
// DOF_Numberer to assign equation numbers to the DOFs.
//
// What: "@(#) DOF_GroupGraph.C, revA"


#include <DOF_GroupGraph.h>
#include <Vertex.h>
#include <AnalysisModel.h>
#include <FE_Element.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>

// constructs the Graph
DOF_GroupGraph::DOF_GroupGraph(AnalysisModel &theModel)
:Graph(theModel.getNumDOF_Groups()+START_VERTEX_NUM), 
 myModel(theModel)
{

    int numVertex = myModel.getNumDOF_Groups();

    if (numVertex == 0) {
	opserr << "WARNING DOF_GroupGraph::DOF_GroupGraph";
	opserr << "  - 0 vertices, has the Domain been populated?\n";
	return;
    }	
	
    DOF_Group *dofPtr;

    // now create the vertices with a reference equal to the DOF_Group number.
    // and a tag which ranges from 0 through numVertex-1

    DOF_GrpIter &dofIter2 = theModel.getDOFs();
    int count = START_VERTEX_NUM;
    while ((dofPtr = dofIter2()) != 0) {
	int DOF_GroupTag = dofPtr->getTag();
	int DOF_GroupNodeTag = dofPtr->getNodeTag();
	int numDOF = dofPtr->getNumFreeDOF();
	Vertex *vertexPtr = new Vertex(DOF_GroupTag, DOF_GroupNodeTag, 0, numDOF);

	if (vertexPtr == 0) {
	    opserr << "WARNING DOF_GroupGraph::DOF_GroupGraph";
	    opserr << " - Not Enough Memory to create ";
	    opserr << count << "th Vertex\n";
	    return;
	}
	
	this->addVertex(vertexPtr);
    }


    // now add the edges, by looping over the Elements, getting their
    // IDs and adding edges between DOFs for equation numbers >= START_EQN_NUM
    
    FE_Element *elePtr;
    FE_EleIter &eleIter = myModel.getFEs();

    while((elePtr = eleIter()) != 0) {
	const ID &id = elePtr->getDOFtags();
	int size = id.Size();
	for (int i=0; i<size; i++) {
	    int dof1 = id(i);
	    for (int j=0; j<size; j++) 
		if (i != j) {
		    int dof2 = id(j);
		    this->addEdge(dof1,dof2);
		}
	}
    }
}

DOF_GroupGraph::~DOF_GroupGraph()
{

}    


