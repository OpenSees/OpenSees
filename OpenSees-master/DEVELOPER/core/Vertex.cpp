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
                                                                        
// $Revision: 1.6 $
// $Date: 2009-05-11 21:37:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/graph/Vertex.cpp,v $
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class implementation for Vertex.
// Vertex is an element of a graph.
//
// What: "@(#) Vertex.C, revA"

#include <Vertex.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>

Vertex::Vertex(int tag, int ref, double weight, int color)
:TaggedObject(tag), myRef(ref), myWeight(weight), myColor(color), 
 myDegree(0), myTmp(0), myAdjacency(0, 8)
{

}    


Vertex::Vertex(const Vertex &other) 
:TaggedObject(other.getTag()), myRef(other.myRef), myWeight(other.myWeight), myColor(other.myColor), 
 myDegree(other.myDegree), myTmp(0), myAdjacency(other.myAdjacency)
{

}

Vertex::~Vertex()
{

}    

void 
Vertex::setWeight(double newWeight) 
{ 
    myWeight = newWeight;
}

void 
Vertex::setColor(int newColor) 
{
    myColor = newColor;
}

void 
Vertex::setTmp(int newTmp) 
{
    myTmp = newTmp;
}

int 
Vertex::getRef(void) const 
{
    return myRef;
}

double
Vertex::getWeight(void) const 
{
    return myWeight;
}

int 
Vertex::getColor(void) const
{
    return myColor;
}

int 
Vertex::getTmp(void) const
{
    return myTmp;
}

int 
Vertex::addEdge(int otherTag)
{
    // don't allow itself to be added
    if (otherTag == this->getTag())
      return 0;

    // check the otherVertex has not already been added
    return myAdjacency.insert(otherTag);
}


int 
Vertex::getDegree(void) const
{
    return myDegree;
}

const ID &
Vertex::getAdjacency(void) const
{
    return myAdjacency;
}

void
Vertex::Print(OPS_Stream &s, int flag)
{
    s << this->getTag() << " " ;
    s << myRef << " ";
    if (flag == 1) 
	s << myWeight << " " ;
    else if (flag ==2) 
	s << myColor << " " ;
    else if (flag == 3)
        s << myWeight << " " << myColor << " " ;    
    else if (flag == 4)
      s << " weight: " << myWeight << " color: " << myColor << " tmp: " << myTmp << " " ;

    s << "ADJACENCY: " << myAdjacency;    	
}


int 
Vertex::sendSelf(int commitTag, Channel &theChannel)
{
  // send the tag/ref/color/degree/tmp, an indication if weighted & size of adjacency
  static ID idData(7);
  idData(0) = this->getTag();
  idData(1) = myRef;
  idData(2) = myColor;
  idData(3) = myDegree;
  idData(4) = myTmp;
  if (myWeight == 0.0)
    idData(5) = 0;
  else
    idData(5) = 1;
  idData(6) = myAdjacency.Size();

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "Graph::sendSelf() - failed to receive the initial data\n";
    return -1;
  }

  // if weighted, send the weight
  if (myWeight != 0.0) {
    static Vector vectData(1);
    vectData(0) = myWeight;
    if (theChannel.sendVector(0, commitTag, vectData) < 0) {
      opserr << "Graph::rendSelf() - failed to receive the weight\n";
      return -2;
    }
  }    

  // finally send the adjacency
  if (theChannel.sendID(0, commitTag, myAdjacency) < 0) {
    opserr << "Graph::sendSelf() - failed to receive the adjacency data\n";
    return -1;
  }  

  return 0;
}


int 
Vertex::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  // recv the tag/ref/color/degree/tmp, an indication if weighted & size of adjacency
  static ID idData(7);
  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "Graph::recvSelf() - failed to receive the initial data\n";
    return -1;
  }
  this->setTag(idData(0));
  myRef = idData(1);
  myColor = idData(2);
  myDegree = idData(3);
  myTmp = idData(4);

  // if weighted, receive the weight
  if (idData(5) == 1) {
    static Vector vectData(1);
    if (theChannel.recvVector(0, commitTag, vectData) < 0) {
      opserr << "Graph::recvSelf() - failed to receive the weight\n";
      return -2;
    }
    myWeight = vectData(0);
  }    

  // resize the adjacency & receive it
  //  myAdjacency[idData(6)-1] = 0;
  int *adjacencyData;
  adjacencyData = new int[idData[6]];
  myAdjacency.setData(adjacencyData, idData[6], true);

  if (theChannel.recvID(0, commitTag, myAdjacency) < 0) {
    opserr << "Graph::recvSelf() - failed to receive the initial data\n";
    return -3;
  }
  return 0;
}
