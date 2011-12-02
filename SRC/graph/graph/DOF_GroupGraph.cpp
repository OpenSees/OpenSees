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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:01:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/graph/DOF_GroupGraph.cpp,v $
                                                                        
                                                                        
// File: ~/graph/graph/DOF_GroupGraph.C
// 
// Written: fmk 
// Created: Sun Sept 15 11:47:47: 1996
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
	
    // create another vertices array which aids in adding edges
    
    int *theDOF_GroupTagVertices = 0;
    int maxDofNum = 0;
    DOF_Group *dofPtr;
    DOF_GrpIter &dofIter = theModel.getDOFs();
    while ((dofPtr = dofIter()) != 0)
	if (dofPtr->getTag() > maxDofNum)
	    maxDofNum = dofPtr->getTag();

    theDOF_GroupTagVertices = new int [maxDofNum+1];

    if (theDOF_GroupTagVertices == 0) {
	opserr << "WARNING DOF_GroupGraph::DOF_GroupGraph ";
	opserr << " - Not Enough Memory for DOF_GroupTagVertices\n";
	return;
    }
    
    for (int j=0; j<=maxDofNum; j++) theDOF_GroupTagVertices[j] = -1;

    // now create the vertices with a reference equal to the DOF_Group number.
    // and a tag which ranges from 0 through numVertex-1

    DOF_GrpIter &dofIter2 = theModel.getDOFs();
    int count = START_VERTEX_NUM;
    while ((dofPtr = dofIter2()) != 0) {
	int DOF_GroupTag = dofPtr->getTag();
	Vertex *vertexPtr = new Vertex(count,DOF_GroupTag);

	if (vertexPtr == 0) {
	    opserr << "WARNING DOF_GroupGraph::DOF_GroupGraph";
	    opserr << " - Not Enough Memory to create ";
	    opserr << count << "th Vertex\n";
	    delete [] theDOF_GroupTagVertices;
	    return;
	}
	
	this->addVertex(vertexPtr);
	theDOF_GroupTagVertices[DOF_GroupTag] = count++;
	
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
	    int vertexTag1 = theDOF_GroupTagVertices[dof1];	    

	    for (int j=0; j<size; j++) 
		if (i != j) {
		    
		    int dof2 = id(j);
		    int vertexTag2 = theDOF_GroupTagVertices[dof2]; 
		    
		    this->addEdge(vertexTag1,vertexTag2);
		}
	}
    }

    // done now delete theDOF_GroupTagVertices
   
    delete [] theDOF_GroupTagVertices;
}

DOF_GroupGraph::~DOF_GroupGraph()
{

}    


