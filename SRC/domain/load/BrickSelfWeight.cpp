///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              BrickSelfWeight
// CLASS:             
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           Finite Element Class
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Zhaohui Yang
// PROGRAMMER:        Zhaohui Yang 
// DATE:              March 2002
// UPDATE HISTORY:
//
//
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef BRICKSELFWEIGHT_CPP
#define BRICKSELFWEIGHT_CPP

                                                                        
// Written: ZHYang UC Davis
// Purpose: This file contains the class definition for 8 node brick self weigth load.

#include <BrickSelfWeight.h>
#include <Vector.h>

Vector BrickSelfWeight::data(1);

BrickSelfWeight::BrickSelfWeight(int tag, const ID &theElementTags)
  :ElementalLoad(tag, LOAD_TAG_BrickSelfWeight, theElementTags)
{

}

BrickSelfWeight::BrickSelfWeight()
  :ElementalLoad(LOAD_TAG_BrickSelfWeight)
{

}

BrickSelfWeight::~BrickSelfWeight()
{

}

const Vector &
BrickSelfWeight::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_BrickSelfWeight;
  //data(0) = P;
  //data(1) = x;
  return data;
}

int 
BrickSelfWeight::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
BrickSelfWeight::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
BrickSelfWeight::Print(OPS_Stream &s, int flag)
{
  s << "BrickSelfWeight...";
  s << "  elements acted on: " << this->getElementTags();
}

#endif

