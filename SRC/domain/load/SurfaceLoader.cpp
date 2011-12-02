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
                                                                        
// Written: Chris McGann, U.Washington
//          02.2011
//
// Description: This file contains the implementation for the SurfaceLoader class.

#include <SurfaceLoader.h>
#include <Vector.h>

Vector SurfaceLoader::data(1);

SurfaceLoader::SurfaceLoader(int tag, int theElementTag)
  : ElementalLoad(tag, LOAD_TAG_SurfaceLoader, theElementTag)
{
}

SurfaceLoader::SurfaceLoader()
  : ElementalLoad(LOAD_TAG_SurfaceLoader)
{
}

SurfaceLoader::~SurfaceLoader()
{
}

const Vector &
SurfaceLoader::getData(int &type, double loadFactor)
{
	type = LOAD_TAG_SurfaceLoader;

	return data;
}

int
SurfaceLoader::sendSelf(int commitTag, Channel &theChannel)
{
	return -1;
}

int
SurfaceLoader::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	return -1;
}

void
SurfaceLoader::Print(OPS_Stream &s, int flag)
{
	s << "SurfaceLoader...";
	s << "  element acted on: " << eleTag << endln;
}

