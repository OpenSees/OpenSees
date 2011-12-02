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
// $Date: 2000-09-15 08:23:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/PlaneFrame.cpp,v $
                                                                        
                                                                        
// File: ~/model/PlaneFrame.C
//
// Written: fmk 
// Created: Mon Sept 15 14:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for PlaneFrame
// PlaneFrame is a class used to model a structure with a plane frame. 
// The object creates the components of the model and adds these to the
// Domain with which it is associated.

#include <PlaneFrame.h>
#include <Domain.h>
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>

#include <ElasticBeam2d.h>
#include <LinearCrdTransf2d.h>
#include <beam2d03.h>
#include <NodalLoad.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <MP_Constraint.h>
#include <LoadPattern.h>
#include <LinearSeries.h>

//  PlaneFrameModel();
//	constructor
PlaneFrame::PlaneFrame(Domain &theDomain)
:ModelBuilder(theDomain)
{

}

// ~PlaneFrame();    
//	destructor,

PlaneFrame::~PlaneFrame()
{

}    

int
PlaneFrame::buildFE_Model(void) 
{
    int numEle, numNode, numSPs, numMPs, numNodLoads;
    int res = 0;
    Domain *theDomain = this->getDomainPtr();    
    if (theDomain == 0) {
	cerr << "FATAL::PlaneFrame::buildModel(void) -";
	cerr << " no associated domain !! - what have you been doing\n";
	exit(-1);
    }	
    
    // get the input file and open it for reading
    char fileName[15];
    cout << "Enter File Containing Model: ";
    cin >> fileName;
    cout << fileName;
    
    ifstream inputFile(fileName, ios::in);
    if (!inputFile) {
	cerr << "FATAL:PlaneFrame::buildModel(void) -";
	cerr << " could not open file: " << fileName << endl;
	exit(-1);
    }
    
    //
    // read in the model paramaters
    // 	  :numNodes numElements numSP_Constraints Num MP_Constraints 
    //     numLoadCases

    inputFile >> numNode >> numEle >> numSPs >> numMPs >> numNodLoads;

    //
    // read in the nodal data
    //  each of numNode nodes: tag xCrd yCrd

    Node *NodePtr;
    int tag, ndof;
    double xCrd, yCrd;
    bool result;
    ndof = 3;
    int i;
	Vector dummy(2);
	CrdTransf2d *theTransf = new LinearCrdTransf2d (1,dummy,dummy);

    
    // for each node read in the data, create a node & add it to the domain
        for (i=0; i<numNode; i++) {
	inputFile >> tag >> xCrd >> yCrd;
	NodePtr = new Node(tag,ndof,xCrd,yCrd);
	result = theDomain->addNode(NodePtr);
        if (result == false) {
                res =-1;
		cerr << "PlaneFrame::buildModel(void) -";
		cerr << " problems adding node " << tag << endl;
        }
    }


    // read in the elemental data, model only recognises type 2 and 3
    // for element type beam2d04 and beam2d03.
    //	  each of numEle elements: type tag A E I node1 node2

    Element *elePtr;
    double A,E,I;
    int nd1, nd2, type;
    for (i=0; i<numEle; i++) {
	inputFile >> type >> tag >> A >> E >> I >> nd1 >> nd2;
	if (type == 2) {
	    elePtr = new ElasticBeam2d(tag,A,E,I,nd1,nd2,*theTransf);
	} else if (type == 3) {
	    elePtr = new beam2d03(tag,A,E,I,nd1,nd2);
	} else {
            res =-1;
	    cerr << "ERROR PlaneFrame::PlaneFrame - Invalid element type: ";
	    cerr << type << endl;
	    return -1;
	}
	result = theDomain->addElement(elePtr);
        if (result == false) {
                res =-1;
		cerr << "PlaneFrame::buildModel(void) -";
		cerr << " problems adding element " << tag << endl;
        }
    }

    // add a LoadPattern with a LinearSeries
    LoadPattern *theLoadPattern = new LoadPattern(0);
    TimeSeries  *theSeries = new LinearSeries(1.0);
    theLoadPattern->setTimeSeries(theSeries);
    result = theDomain->addLoadPattern(theLoadPattern);
    if (result == false) {
      cerr << "PlaneFrame::buildModel(void) -";
      cerr << " problems adding load pattern " << *theLoadPattern;
      res =-1;
    }    

    
    //
    // read in the SP Constraint data
    //   each of numSPs SP_Constraints: nodeTag dof value
    
    SP_Constraint *SPPtr;
    for (i=0; i<numSPs; i++) {
	inputFile >> nd1 >> nd2 >> E;
	SPPtr = new SP_Constraint(i,nd1,nd2,E);
	result = theDomain->addSP_Constraint(SPPtr, 0);
        if (result == false) {
                res =-1;
		cerr << "PlaneFrame::buildModel(void) -";
		cerr << " problems adding SP_Constraint on " << nd1 << endl;
        }
    }    

    //
    // read in the MP Constraint data
    //   each of numMPs SP_Constraints: nodeTag dof value
    
    MP_Constraint *MPPtr;
    int numDOF1, numDOF2;
    for (i=0; i<numMPs; i++) {
	inputFile >> nd1 >> nd2 >> numDOF1 >> numDOF2;
	ID constrainedDOF(numDOF1);
	ID retainedDOF(numDOF2);	
	inputFile >> constrainedDOF;
	inputFile >> retainedDOF;
	Matrix Ccr(numDOF1,numDOF2);
	inputFile >> Ccr;
	MPPtr = new MP_Constraint(i,nd2,nd1,Ccr,constrainedDOF,retainedDOF);
	result = theDomain->addMP_Constraint(MPPtr);
        if (result == false) {
                res =-1;
		cerr << "PlaneFrame::buildModel(void) -";
		cerr << " problems adding MP_Constraint on " << nd1 << endl;
        }

    }        

    Vector forces(3);
    NodalLoad *nodeLoadPtr;    
    for (i=0; i<numNodLoads; i++) {
	inputFile >> tag >> forces(0) >> forces(1) >> forces(2);
	nodeLoadPtr = new NodalLoad(i, tag, forces);
	bool result = theDomain->addNodalLoad(nodeLoadPtr, 0);
	if (result == false) {
		cerr << "PlaneFrame::buildModel(void) -";
		cerr << " problems adding load " << *nodeLoadPtr;
                res =-1;
        }
    }    

    //
    // done inputing the model
    //
    return res;
}	
	





