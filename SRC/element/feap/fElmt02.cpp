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
// $Date: 2000-09-15 08:23:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/feap/fElmt02.cpp,v $
                                                                        
                                                                        
// File: ~/element/fortran/fElmt02.C
// 
// Written: fmk 
// Created: 03/99
// Revision: A
//
// Description: This file contains the implementation for the fElmt02 class.
//
// What: "@(#) fElement.C, revA"

#include "fElmt02.h"
#include <ID.h>
#include <Vector.h>
#include <elementAPI.h>
#include <Channel.h>

void* OPS_fElmt02()
{

    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();

    if (ndm != 2 && ndf != 2) {
	opserr << "WARNING - fTruss eleTag? iNode? jNode? A? E? needs ndm=2, ndf=2\n";
	return 0;
    }

    // check the number of arguments is correct
    int argc = OPS_GetNumRemainingInputArgs() + 2;
    int eleArgStart = 2;
    if ((argc-eleArgStart) < 5) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: element fTruss eleTag? iNode? jNode? A? E?\n";
	return 0;
    }    


    // get the id and end nodes
    int idata[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid truss eleTag, iNode or jNode" << endln;
	return 0;
    }
    int trussId=idata[0], iNode=idata[1], jNode=idata[2];

    double ddata[2];
    numdata = 2;
    if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
	opserr << "WARNING invalid truss A or E" << endln;
	return 0;
    }
    
    double A=ddata[0],E=ddata[1];
 
    // now create the truss and add it to the Domain
    return new fElmt02(trussId,iNode,jNode,A,E);
}

fElmt02::fElmt02(int tag, int nd1, int nd2, double A, double E, double rho)
:fElement(tag, ELE_TAG_fElmt02, 2, 3, 2, 2, 2, 0, 0)
{
    (*data)(0) = A;
    (*data)(1) = E;
    (*data)(2) = rho;
    
    (*connectedNodes)(0) = nd1; 
    (*connectedNodes)(1) = nd2;   
}

fElmt02::fElmt02(int tag, int nd1, int nd2, int iow)
:fElement(tag, ELE_TAG_fElmt02, 2, 3, 2, 2, 2, iow)
{
    (*connectedNodes)(0) = nd1; 
    (*connectedNodes)(1) = nd2;   
}
    
fElmt02::fElmt02()
:fElement(ELE_TAG_fElmt02)    
{
    // does nothing
}

fElmt02::~fElmt02()
{
    // does nothing
}

int 
fElmt02::sendSelf(int commitTag, Channel &theChannel)
{
  return this->fElement::sendSelf(commitTag, theChannel);
}
int 
fElmt02::recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  return this->fElement::recvSelf(commitTag, theChannel, theBroker);
}
