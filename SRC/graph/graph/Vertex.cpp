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
// $Source: /usr/local/cvs/OpenSees/SRC/graph/graph/Vertex.cpp,v $
                                                                        
                                                                        
// File: ~/graph/graph/Vertex.C
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class implementation for Vertex.
// Vertex is an element of a graph.
//
// What: "@(#) Vertex.C, revA"

#include <Vertex.h>

Vertex::Vertex(int tag, int ref, double weight, int color)
:TaggedObject(tag), myRef(ref), myWeight(weight), myColor(color), 
 myDegree(0), myTmp(0), myAdjacency(0,8)
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
    if (myAdjacency.getLocation(otherTag) < 0) {
	myAdjacency[myDegree]  = otherTag;
	myDegree++;
	return 0;
    }
    else
	return 1;
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
Vertex::Print(ostream &s, int flag)
{
    s << this->getTag() << " " ;
    s << myRef << " ";
    if (flag == 1) 
	s << myWeight << " " ;
    else if (flag ==2) 
	s << myColor << " " ;
    else if (flag == 3)
	s << myWeight << " " << myColor << " " ;    

    s << "ADJACENCY: " << myAdjacency;    	
}


