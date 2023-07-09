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
                                                                        
// $Revision: 1.5 $
// $Date: 2005-11-14 20:47:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/graph/DOF_Graph.cpp,v $
                                                                        
                                                                        
// File: ~/graph/graph/DOF_Graph.C
// 
// Written: fmk 
// Created: Sun Sept 15 11:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for DOF_Graph.
// DOF_Graph is a graph of the DOFs in the analysis model. It is used
// by the SysOfEqn to determine its size.
//
// What: "@(#) DOF_Graph.C, revA"


#include <DOF_Graph.h>
#include <Vertex.h>
#include <AnalysisModel.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <FE_Element.h>
#include <FE_EleIter.h>

#define START_EQN_NUM 0
// constructs the Graph
// assumes eqn numbers are numbered continuously from START_EQN_NUM

DOF_Graph::DOF_Graph(AnalysisModel &theModel)
:Graph(theModel.getNumEqn()), 
 myModel(theModel)
{

  /********************** old way to create vertices with dof Tags ********
    int numVertex = myModel.getNumEqn();

    if (numVertex <= 0) {
	opserr << "WARNING DOF_Graph::DOF_Graph";
	opserr << "  - 0 equations?\n";
	return;
    }	

    // now create the vertices with a reference equal to the eqn number.
    // and a tag which ranges from START_VERTEX_NUM through 

    for (int i =0; i<numVertex; i++) {
	Vertex *vertexPtr = new Vertex(i, i);
	
	if (vertexPtr == 0) {
	    opserr << "WARNING DOF_Graph::DOF_Graph";
	    opserr << " - Not Enough Memory to create " << i+1 << "th Vertex\n";
	    return;
	}
	this->addVertex(vertexPtr,false);	
    }
    *****************************************************************************/

  //
  // create a vertex for each dof
  //

  DOF_Group *dofPtr =0;
  DOF_GrpIter &theDOFs = myModel.getDOFs();
  while ((dofPtr = theDOFs()) != 0) {
    const ID &id = dofPtr->getID();
    int size = id.Size();
    for (int i=0; i<size; i++) {
      int dofTag = id(i);
      if (dofTag >= START_EQN_NUM) {
	Vertex *vertexPtr = this->getVertexPtr(dofTag);
	if (vertexPtr == 0) {
	  Vertex *vertexPtr = new Vertex(dofTag, dofTag);      
	  if (vertexPtr == 0) {
	    opserr << "WARNING DOF_Graph::DOF_Graph";
	    opserr << " - Not Enough Memory to create " << i+1 << "th Vertex\n";
	    return;
	  }
	  if (this->addVertex(vertexPtr, false) == false) {
	    opserr << "WARNING DOF_Graph::DOF_Graph - error adding vertex\n";
	  }
	}
      }
    }
  }

  // now add the edges, by looping over the FE_elements, getting their
  // IDs and adding edges between DOFs for equation numbers >= START_EQN_NUM
  
  FE_Element *elePtr =0;
  FE_EleIter &eleIter = myModel.getFEs();
  int cnt = 0;

  while((elePtr = eleIter()) != 0) {
    const ID &id = elePtr->getID();
    cnt++;
    int size = id.Size();
    for (int i=0; i<size; i++) {
      int eqn1 = id(i);
	    
      // if eqnNum of DOF is a valid eqn number add an edge
      // to all other DOFs with valid eqn numbers.
      
      if (eqn1 >=START_EQN_NUM) {
	for (int j=i+1; j<size; j++) {
	  int eqn2 = id(j);
	  if (eqn2 >=START_EQN_NUM)
	    this->addEdge(eqn1-START_EQN_NUM+START_VERTEX_NUM,
			  eqn2-START_EQN_NUM+START_VERTEX_NUM);
	}
      }
    }
  }
}

DOF_Graph::~DOF_Graph()
{

}    




