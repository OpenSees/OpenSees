/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** See file 'COPYRIGHT'  in main directory for information on usage   **
** and redistribution of OpenSees, and for a DISCLAIMER OF ALL        **
** WARRANTIES.                                                        **
**                                                                    **
** Element2dGNL.cpp: implementation of the Element2dGNL class         **
** Developed by:                                                      **
**    Rohit Kaul       (rkaul@stanford.edu)                           **
**    Greg Deierlein   (ggd@stanford.edu)                             **
**                                                                    **
**           John A. Blume Earthquake Engineering Center              **
**                    Stanford University                             **
** ****************************************************************** **/

#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <math.h>
#include <stdlib.h>
#include <elementAPI.h>

#include "Elastic2DGNL.h"

#define Ele_TAG_Elastic2dGNL -1
 
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void * OPS_ADD_RUNTIME_VPV(OPS_Elastic2DGNL)
{
    if (OPS_GetNumRemainingInputArgs() < 6)
    {
	opserr << "WARNING insufficient arguments\n";
	opserr << "element element2dGNL int tag, int Nd1, int Nd2, double A, double E, double Iz, <int linear>\n";

	return 0;
    }

    int idata[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid Elastic2dGNL int inputs" << endln;
	return 0;
    }
    
    int tag = idata[0];
    int ndI = idata[1];
    int ndJ = idata[2];

    double data[3];
    numdata = 3;
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid Elastic2dGNL double inputs" << endln;
	return 0;
    }

    double A = data[0];
    double E = data[1];
    double I = data[2];
    
    bool   linear = false;

    if (OPS_GetNumRemainingInputArgs() > 0) {
	numdata = 1;
	if (OPS_GetIntInput(&numdata, idata) < 0) {
	    opserr << "WARNING invalid Elastic2dGNL int inputs" << endln;
	    return 0;
	}
	if (idata[0] == 1) linear = true;
    }
	
    return new Elastic2dGNL(tag, A, E, I, ndI, ndJ, linear);//, false, massDens);
}

// Constructor
Elastic2dGNL::Elastic2dGNL(int tag, double a, double e, double i, int Nd1, int Nd2, 
			   bool islinear, double rho)
  :UpdatedLagrangianBeam2D(tag, Ele_TAG_Elastic2dGNL, Nd1, Nd2, islinear),
   A(a), E(e), Iz(i)
{
  massDof = A*L*rho;
  massDof = massDof/2;

}

Elastic2dGNL::~Elastic2dGNL()
{

}

void Elastic2dGNL::getLocalMass(Matrix &M)
{
    if(massDof < 0)
    {
            opserr << "Elastic2dGNL::getMass - Distributed mass not implemented\n";
	    M.Zero();
    }
    else if(massDof == 0)//this cond. is taken care of already
    {
        M.Zero();
    }
    else
    {
        M.Zero();
	    M(0,0) = M(1,1) = M(2,2) = M(3,3) = M(4,4) = M(5,5) = massDof;
    }

}

void Elastic2dGNL::getLocalStiff(Matrix &K)
{
 double	EIbyL = E*Iz/L_hist;
 double l = L_hist;

    K(0, 1) = K(0, 2) = K(0, 4) = K(0, 5)=0;
    K(1, 0) = K(1, 3) =0;
    K(2, 0) = K(2, 3) =0;
    K(3, 1) = K(3, 2) = K(3, 4) = K(3, 5)=0;
    K(4, 0) = K(4, 3) =0;
    K(5, 0) = K(5, 3) =0;

	K(0,0) = K(3,3) = (A/Iz)*(EIbyL);
	K(0,3) = K(3,0) = (-A/Iz)*(EIbyL);
	K(1,1) = K(4,4) = (12/(l*l))*(EIbyL);
	K(1,4) = K(4,1) = (-12/(l*l))*(EIbyL);
	K(1,2) = K(2,1) = K(1,5) = K(5,1) = (6/l)*(EIbyL);
	K(2,4) = K(4,2) = K(4,5) = K(5,4) = (-6/l)*(EIbyL);
	K(2,2) = K(5,5) = 4*(EIbyL);
	K(2,5) = K(5,2) = 2*(EIbyL);
	

}//getLocalStiff



void Elastic2dGNL::Print(OPS_Stream &s, int flag)
{
    s << "\nElement No: " << this->getTag();
    s << " type: Elastic2dGNL  iNode: " << connectedExternalNodes(0);
    s << " jNode: " << connectedExternalNodes(1);
	if(isLinear) s << "(1st-Order):\n";
	else		 s << "(2nd-Order):\n";
}

int Elastic2dGNL::sendSelf(int commitTag, Channel &theChannel)
{
	opserr << "WARNING (W_C_10) - Elastic2dGNL::sendSelf(..) [" << getTag() <<"]\n";
	opserr << "method not implemented\n";
	return -1;
}

int Elastic2dGNL::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
        opserr << "WARNING (W_C_20) - Elastic2dGNL::recvSelf(..) [" << getTag() <<"]\n";
        opserr << "method not implemented\n";
	return -1;
}

