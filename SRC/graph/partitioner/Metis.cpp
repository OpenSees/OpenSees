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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/partitioner/Metis.cpp,v $
                                                                        
                                                                        
// File: ~/graph/partitioner/Metis.C
// 
// Written: fmk 
// Created: Sun Sept 15 11:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for Metis.
// Metis is a type of GraphPartitioner which uses 'METIS - Unstructured
// Graph Partitioning And Sparse Matrix Ordering System', developed by
// G. Karypis and V. Kumar at the University of Minnesota. The metis
// files are found in metis-2.0 which were downloaded.
//     This class provides the C++ interface for metis which will allow
// it to fit seamlessly into our system.
//
// What: "@(#) Metis.C, revA"

#include <Metis.h>
#include <Graph.h>
#include <Vertex.h>

/* stuff needed to get the program working on the clump & NOW machines*/
#include <bool.h>
#ifdef _LINUX
extern "C" {
#include <GKlib.h>
}

int IsWeighted;
timer TotalTmr;         /* Times the entire algorithm */
timer CoarsenTmr;       /* Times total coarsening time */
timer GreedyTmr;        /* Times Total Greedy Time */
timer GreedyInitTmr;    /* Times initialization cost of Greedy */
timer GreedyIterTmr;    /* Times Iterative cost of Greedy */
timer GreedyWrapUpTmr;  /* Times the wrap up phase of Greedy */
timer MlevelTmr;        /* Times the entire multilevel algorithm */
timer InitPartTmr;      /* Times the initial partition phase */
timer ProjectTmr;       /* Times the projection of the partition */
timer SplitTmr;         /* Times the boundary creation */
timer BalanceTmr;       /* Times the time required for balancing */
timer IOTmr;            /* Times the file input time */
timer UncrsTmr;         /* Times the file input time */
#endif

extern "C" 
int PMETIS(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);

extern "C" 
int KMETIS(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);


Metis::Metis(int numParts) 
:GraphNumberer(GraphNUMBERER_TAG_Metis),
 myPtype(0), myMtype(0), myCoarsenTo(0), myRtype(0), myIPtype(0),
 defaultOptions(true), numPartitions(numParts), theRefResult(0)
{
    
}

Metis::Metis(int Ptype, 
	  int Mtype, 
	  int coarsenTo,
	  int Rtype, 
      	  int IPtype,
	  int numParts)
:GraphNumberer(GraphNUMBERER_TAG_Metis),
 myPtype(Ptype), myMtype(Mtype), myCoarsenTo(coarsenTo), myRtype(Rtype), 
 myIPtype(IPtype), defaultOptions(false), numPartitions(numParts), theRefResult(0)
{
    // check the options are valid
    checkOptions();
}    

    
Metis::~Metis() 
{

}    

bool
Metis::setOptions(int Ptype, 
		  int Mtype,
		  int coarsenTo,
		  int Rtype, 
		  int IPtype)
{
    myPtype=Ptype; 
    myMtype = Mtype; 
    myCoarsenTo = coarsenTo; 
    myRtype = Rtype; 
    myIPtype = IPtype; 

    defaultOptions = false;
    
    return checkOptions();
    
}    


// bool checkOptions(void) const
//	returns true if options are o.k., false otherwise

bool
Metis::checkOptions(void) 
{

    // check default not set
    
    if (defaultOptions == true)
	return true;

    
    // otherwise check all options for valid value
    bool okFlag = true;    
    
    if ((myPtype != 1) || (myPtype != 2)) {
	okFlag = false;
	cerr << "WARNING: Metis::partition ";
	cerr << " - Illegal Ptype " << myPtype << endl;
    }
	
    if ((myMtype != 1) ||  (myMtype != 2) || (myMtype != 3) || 
	((myMtype != 4) || (myMtype != 5) || myMtype != 11) || 
	(myMtype != 21) || (myMtype != 51)) {
	okFlag = false;
	cerr << "WARNING: Metis::partition ";
	cerr << " - Illegal Mtype " << myMtype << endl;
    }

    if (myPtype == 1)  
	if ((myRtype != 1) || (myRtype != 2) || (myRtype != 3) ||
	    (myRtype != 11) || (myRtype != 12) || (myRtype != 13) ||
	    (myRtype != 20)) {
	    okFlag = false;
	    cerr << "WARNING: Metis::partition ";
	    cerr << " - Illegal Rtype " << myRtype << endl;
	    cerr << " for Ptype " << myPtype << endl;	    
	}
    else 
	if (myPtype == 2)
	    if ((myRtype != 11) || (myRtype != 12) || (myRtype != 20)) {
		okFlag = false;
		cerr << "WARNING: Metis::partition ";
		cerr << " - Illegal Rtype " << myRtype << endl;
		cerr << " for Ptype " << myPtype << endl;	    
	    }
    
    if ((myIPtype != 1) || (myIPtype != 2) || (myIPtype != 3) || 
	(myIPtype != 4)) {
	okFlag = false;
	cerr << "WARNING: Metis::partition ";
	cerr << " - Illegal IPtype " << myIPtype << endl;
    }	    
	
    if (myCoarsenTo < 0) {
	okFlag = false;
	cerr << "WARNING: Metis::partition ";
	cerr << " - Illegal coarsen To " << myCoarsenTo << endl;
    }	    

    if (okFlag == false)
	defaultOptions = true;
    
    return okFlag;

}


bool 
Metis::setDefaultOptions(void)
{
    defaultOptions = true;
    return true;
}


// int partition(Graph &theGraph, int numPart)
//	Method to partition the graph. It first creates the arrays needed
//	by the metis lib and then invokes a function from the metis lib to
//	partition the graph. The solors of the vertices of the graph are
//	set to colors 0 through numPart-1 to indicate which partition the
//	vrtices are in. Returns -1 if options are not set, -2 if metis failed.

int
Metis::partition(Graph &theGraph, int numPart)
{
    // first we check that the options are valid
    if (checkOptions() == false)
	return -1;
    
    // now we get room for the data structures metis needs

    int numVertex = theGraph.getNumVertex();
    int numEdge = theGraph.getNumEdge();

    int *options = new int [5];
    int *partition = new int [numVertex+1];    
    int *xadj = new int [numVertex+2];
    int *adjncy = new int [2*numEdge];
    int *vwgts = 0;
    int *ewgts = 0;
    int numbering = 0;
    int weightflag = 0; // no weights on our graphs yet

    if (START_VERTEX_NUM == 0)
	numbering = 0;	
    else if (START_VERTEX_NUM == 1)
	numbering = 1;
    else {
	cerr << "WARNING Metis::partition - No partitioning done";
	cerr << " vertex numbering must start at 0 or 1\n";
	return (-2);
    }
    int edgecut;
    
    if ((options == 0) || (partition == 0) || (xadj == 0) || (adjncy == 0)) {
	cerr << "WARNING Metis::partition - No partitioning done";
	cerr << " as ran out of memory\n";
	return (-2);
    }


    // we build these data structures
    
    int indexEdge = 0;
    xadj[0] = 0;

    Vertex *vertexPtr;
    for (int vertex =0; vertex<numVertex; vertex++) {
	vertexPtr = theGraph.getVertexPtr(vertex+START_VERTEX_NUM);
	
	// check we don't have an invalid vertex numbering scheme
	// if so WARNING message, clean up and return -2

	if (vertexPtr == 0) {
	    cerr << "WARNING Metis::partition - No partitioning done";
	    cerr << " Metis requires consequtive Vertex Numbering\n";
	    
	    delete [] options;
	    delete [] partition;
	    delete [] xadj;
	    delete [] adjncy;
	    
	    return -2;
	}
	
	const ID&adjacency = vertexPtr->getAdjacency();
	int degree = adjacency.Size();
	for (int i=0; i<degree; i++) {
	    adjncy[indexEdge++] = adjacency(i)-START_VERTEX_NUM;
	}
	
	xadj[vertex+1] = indexEdge;
    }


    if (defaultOptions == true) 
	options[0] = 0;
    else {
	options[0] =1;
	options[1] = myCoarsenTo;
	options[2] = myMtype;
	options[3] = myIPtype;
	options[4] = myRtype;
    }
    
    
    // we now the metis routines

    if (myPtype == 1)
	PMETIS(&numVertex, xadj, adjncy, vwgts, ewgts, &weightflag, &numPart,
	       options, &numbering, &edgecut, partition);
    else		
	KMETIS(&numVertex, xadj, adjncy, vwgts, ewgts, &weightflag, &numPart,
	       options, &numbering, &edgecut, partition);


    // we set the vertex colors to correspond to the partitioned scheme
    for (int vert =0; vert<numVertex; vert++) {
	vertexPtr = theGraph.getVertexPtr(vert+START_VERTEX_NUM);
	vertexPtr->setColor(partition[vert]+1); // start colors at 1
    }

    // clean up the space and return
    
    delete [] options;
    delete [] partition;
    delete [] xadj;
    delete [] adjncy;
    
    return 0;
}






const ID &
Metis::number(Graph &theGraph, int lastVertex)
{
    // first we check that the options are valid
    // first check our size, if not same make new
    int numVertex = theGraph.getNumVertex();
    // delete the old
    if (theRefResult != 0)
      delete theRefResult;
	
    theRefResult = new ID(numVertex);

    if (theRefResult == 0) {
      cerr << "ERROR:  Metis::number - Out of Memory\n";
      theRefResult = new ID(0);
      return *theRefResult;
    }

    if (checkOptions() == false) {
      cerr << "ERROR:  Metis::number - chek options failed\n";
      return *theRefResult;
    }
    
    // now we get room for the data structures metis needs
    int numEdge = theGraph.getNumEdge();

    int *options = new int [5];
    int *partition = new int [numVertex+1];    
    int *xadj = new int [numVertex+2];
    int *adjncy = new int [2*numEdge];
    int *vwgts = 0;
    int *ewgts = 0;
    int numbering = 0;
    int weightflag = 0; // no weights on our graphs yet

    if (START_VERTEX_NUM == 0)
	numbering = 0;	
    else if (START_VERTEX_NUM == 1)
	numbering = 1;
    else {
	cerr << "WARNING Metis::partition - No partitioning done";
	cerr << " vertex numbering must start at 0 or 1\n";
	return *theRefResult;
    }
    int edgecut;
    
    if ((options == 0) || (partition == 0) || (xadj == 0) || (adjncy == 0)) {
	cerr << "WARNING Metis::partition - No partitioning done";
	cerr << " as ran out of memory\n";
	return *theRefResult;
    }


    // we build these data structures
    
    int indexEdge = 0;
    xadj[0] = 0;

    Vertex *vertexPtr;
    for (int vertex =0; vertex<numVertex; vertex++) {
	vertexPtr = theGraph.getVertexPtr(vertex+START_VERTEX_NUM);
	
	// check we don't have an invalid vertex numbering scheme
	// if so WARNING message, clean up and return -2

	if (vertexPtr == 0) {
	    cerr << "WARNING Metis::partition - No partitioning done";
	    cerr << " Metis requires consequtive Vertex Numbering\n";
	    
	    delete [] options;
	    delete [] partition;
	    delete [] xadj;
	    delete [] adjncy;
	    
	    return *theRefResult;
	}
	
	const ID&adjacency = vertexPtr->getAdjacency();
	int degree = adjacency.Size();
	for (int i=0; i<degree; i++) {
	    adjncy[indexEdge++] = adjacency(i)-START_VERTEX_NUM;
	}
	
	xadj[vertex+1] = indexEdge;
    }


    if (defaultOptions == true) 
	options[0] = 0;
    else {
	options[0] =1;
	options[1] = myCoarsenTo;
	options[2] = myMtype;
	options[3] = myIPtype;
	options[4] = myRtype;
    }
    
    
    // we now the metis routines

    if (myPtype == 1)

      PMETIS(&numVertex, xadj, adjncy, vwgts, ewgts, &weightflag, 
	     &numPartitions, options, &numbering, &edgecut, partition);
    else		
      KMETIS(&numVertex, xadj, adjncy, vwgts, ewgts, &weightflag, 
	     &numPartitions, options, &numbering, &edgecut, partition);

cerr << "Metis::number -2\n";
    // we assign numbers now based on the partitions returned.
    // each vertex in partion i is assigned a number less than
    // thos in partion i+1: NOTE WE DON'T CARE WHAT THE NUMBERING IS
    // WITHIN A PARTITION
    int count = 0;
    for (int i=0; i<numPartitions; i++) {
      for (int vert=0; vert<numVertex; vert++) {
	if (partition[vert] == i) {
	  vertexPtr = theGraph.getVertexPtr(vert+START_VERTEX_NUM);
	  vertexPtr->setTmp(count+1);
	  (*theRefResult)(count) = vertexPtr->getRef();
	  count++;
	}
      }
    }
cerr << "Metis::number -3\n";
    // clean up the space and return
    delete [] options;
    delete [] partition;
    delete [] xadj;
    delete [] adjncy;
cerr << "Metis::number -4\n";    
    return *theRefResult;
}


const ID &
Metis::number(Graph &theGraph, const ID &lastVertices)
{
  return this->number(theGraph);
}

int 
Metis::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int 
Metis::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}

