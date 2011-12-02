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
// $Date: 2009-05-14 23:25:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/ArpackSOE.cpp,v $

// Written: fmk
// Created: 05/09
//
// Description: This file contains the class definition for ArpackSOE

#include <ArpackSOE.h>
#include <ArpackSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>

ArpackSOE::ArpackSOE(LinearSOE &theLinSOE,
		     double theShift)
:EigenSOE(EigenSOE_TAGS_ArpackSOE),
 M(0), Msize(0), mDiagonal(false), shift(theShift), theModel(0), theSOE(&theLinSOE),
 processID(-1), numChannels(0), theChannels(0), localCol(0), sizeLocal(0)
{
  ArpackSolver *theSolvr = new ArpackSolver();
  this->setSolver(*theSolvr);
  theSolvr->setEigenSOE(*this);
}


ArpackSOE::ArpackSOE()
:EigenSOE(EigenSOE_TAGS_ArpackSOE),
 M(0), Msize(0), mDiagonal(false), shift(0), theModel(0), theSOE(0),
 processID(-1), numChannels(0), theChannels(0), localCol(0), sizeLocal(0)
{
  ArpackSolver *theSolvr = new ArpackSolver();
  this->setSolver(*theSolvr);
  theSolvr->setEigenSOE(*this);
}


int
ArpackSOE::getNumEqn(void) const
{
  return theSOE->getNumEqn();
}
    
ArpackSOE::~ArpackSOE()
{
  if (M != 0) delete [] M;

  if (theSOE != 0)
    delete theSOE;
}

int 
ArpackSOE::setSize(Graph &theGraph)
{
  int result = 0;
  int size = 0;

  if (processID == -1) {

    size = theGraph.getNumVertex();    

  } else {

    // fist itearte through the vertices of the graph to get n
    int maxVertexTag = -1;
    Vertex *theVertex;
    
    VertexIter &theVertices = theGraph.getVertices();
    while ((theVertex = theVertices()) != 0) {
      int vertexTag = theVertex->getTag();
      if (vertexTag > maxVertexTag)
	maxVertexTag = vertexTag;
    }

    if (processID != 0) {
      
      //
      // if subprocess, send local max vertexTag (n)
      // recv ax n from P0
      //

      static ID data(1);
      
      data(0) = maxVertexTag;
      Channel *theChannel = theChannels[0];
      theChannel->sendID(0, 0, data);    // send local max tag
      theChannel->recvID(0, 0, data);    // recv global max tag 

      size = data(0);    
    } else {

      //
      // from each distributed soe recv it's max n and compare; return max n to all
      //
      
      static ID data(1);
      FEM_ObjectBroker theBroker;

      for (int j=0; j<numChannels; j++) {
	Channel *theChannel = theChannels[j];
	theChannel->recvID(0, 0, data);
	if (data(0) > maxVertexTag)
	  maxVertexTag = data(0);
      }
      
      data(0) = maxVertexTag;

      for (int j=0; j<numChannels; j++) {
	Channel *theChannel = theChannels[j];
	theChannel->sendID(0, 0, data);
      }
      size = maxVertexTag;
    }

    size+=1; // vertices numbered 0 through n-1
  }
  
  // invoke setSize on the linearSOE

  result = theSOE->setSize(theGraph);
  if (result < 0) {
    opserr << "WARNING ArpackSOE::ArpackSOE : - LinearSOE - failed in setSize\n";
    return -1;
  }

  if (size != Msize && size > 0) {

    if (M != 0) 
      delete [] M;
    
    M = new double[size];
    
    if (M == 0) {
      opserr << "WARNING ArpackSOE::ArpackSOE : - out of memory creating memory for M\n";
      Msize = 0;
    } else
      Msize = size;
  }

  //
  // invoke setSize() on the Solver
  //

  EigenSolver *theSolvr = this->getSolver();

  if (theSolvr == 0) {
    opserr << "ArpackSOE::setSize(Graph &theGraph) - no EigenSolver set\n";             
    return -1;
  }
  int solverOK = theSolvr->setSize();

  if (solverOK < 0) {
    opserr << "WARNING:ArpackSOE::setSize() -  solver failed setSize()\n";
    return solverOK;
  } 

  return result;    
}

int 
ArpackSOE::addA(const Matrix &m, const ID &id, double fact)
{

  // check for a quick return 
  if (fact == 0.0)  return 0;

  return theSOE->addA(m, id, fact);
}


void 
ArpackSOE::zeroA(void)
{
  return theSOE->zeroA();
}

int 
ArpackSOE::addM(const Matrix &m, const ID &id, double fact)
{
  int res = this->addA(m, id, -shift);

  if (res < 0)
    return res;

  if (mDiagonal == false)
    return  res;

  int idSize = id.Size();
  for (int i=0; i<idSize; i++) {
    int locI = id(i);
    if (locI >= 0 && locI < Msize) {
      for (int j=0; j<idSize; j++) {
	int locJ = id(j);
	if (locJ >= 0 && locJ < Msize) {
	  if (locI == locJ) {
	    M[locI] += m(i,i);
	  } else {
	    if (m(i,j) != 0.0) {
	      mDiagonal = false;
	      return res;
	    }
	  }
	}
      }
    }
  }

  return 0;
}   
 
void 
ArpackSOE::zeroM(void)
{
  mDiagonal = true;

  for (int i=0; i<Msize; i++)
    M[i] = 0;
}


double 
ArpackSOE::getShift(void)
{
    return shift;
}


int 
ArpackSOE::sendSelf(int commitTag, Channel &theChannel)
{
  int sendID =0;
  
  if (processID == -1)
    processID = 0;
  
  // if P0 check if already sent. If already sent use old processID; if not allocate a new process 
  // id for remote part of object, enlarge channel * to hold a channel * for this remote object.
  // if not P0, send current processID

  if (processID == 0) {
    
    // check if already using this object
    bool found = false;
    for (int i=0; i<numChannels; i++)
      if (theChannels[i] == &theChannel) {
	sendID = i+1;
	found = true;
      }
    
    // if new object, enlarge Channel pointers to hold new channel * & allocate new ID
    if (found == false) {
      int nextNumChannels = numChannels + 1;
      Channel **nextChannels = new Channel *[nextNumChannels];
      if (nextNumChannels == 0) {
	opserr << "ArpackSOE::sendSelf() - failed to allocate channel array of size: " << 
	  nextNumChannels << endln;
	return -1;
      }
      for (int i=0; i<numChannels; i++)
	nextChannels[i] = theChannels[i];
      nextChannels[numChannels] = &theChannel;
      
      numChannels = nextNumChannels;
      
      if (theChannels != 0)
	delete [] theChannels;
      
      theChannels = nextChannels;
      
      if (localCol != 0)
	delete [] localCol;
      localCol = new ID *[numChannels];
      if (localCol == 0) {
	opserr << "ArpackSOE::sendSelf() - failed to allocate id array of size: " << 
	  nextNumChannels << endln;
	return -1;
      }
      for (int i=0; i<numChannels; i++)
	localCol[i] = 0;    
      
      if (sizeLocal != 0)
	delete sizeLocal;
      
      sizeLocal = new ID(numChannels);
      
      // allocate new processID for remote object
      sendID = numChannels;
    }
  } else 
    sendID = processID;
  
  // send remotes processID
  ID idData(2);
  idData(0) = sendID;
  idData(1) = theSOE->getClassTag();
  
  int res = theChannel.sendID(0, commitTag, idData);
  if (res < 0) {
    opserr <<"WARNING ArpackSOE::sendSelf() - failed to send data\n";
    return -1;
  }

  res = theSOE->sendSelf(commitTag, theChannel);

  if (res < 0) {
    opserr <<"WARNING ArpackSOE::sendSelf() - failed to send LinearSOE\n";
    return -1;
  }

  return 0;  
}

    
int 
ArpackSOE::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  ID idData(2);
  int res = theChannel.recvID(0, commitTag, idData);
  if (res < 0) {
    opserr <<"WARNING ArpackSOE::recvSelf() - failed to send data\n";
    return -1;
  }	      
  processID = idData(0);

  numChannels = 1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;

  localCol = new ID *[numChannels];
  for (int i=0; i<numChannels; i++)
    localCol[i] = 0;

  if (sizeLocal != 0)
    delete sizeLocal;

  sizeLocal = new ID(numChannels);

  theSOE = theBroker.getNewLinearSOE(idData(1)); 
  if (theSOE == 0) {
    opserr <<"WARNING ArpackSOE::recvSelf() - failed to get LinearSOE of type: " << idData(1) << endln;
    return -1;
  }	      

  if (theSOE->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr <<"WARNING ArpackSOE::recvSelf() - LinearSOE of type: " << idData(1) << " failed in recvSelf\n";
    return -1;
  }

  return 0;
}

int 
ArpackSOE::setLinks(AnalysisModel &theAnalysisModel)
{
  theModel = &theAnalysisModel;
  return 0;
}

int 
ArpackSOE::setLinearSOE(LinearSOE &theLinearSOE)
{
  theSOE = &theLinearSOE;
	return 0;
}
